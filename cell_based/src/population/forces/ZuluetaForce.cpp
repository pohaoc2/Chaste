/*

Copyright (c) 2005-2024, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "ZuluetaForce.hpp"

template<unsigned DIM>
Zuluetaforce<DIM>::Zuluetaforce()
   : AbstractForce<DIM>(),
     mAreaElasticityParameter(0.001), // These parameters are Case I in Farhadifar's paper
     mTargetAreaParameter(1.0),
     mbeta(0.602),
     mmu(0.02),
     mgamma(0.1),
     mmc(0.13),
     mAmax(0.8),
     mK(0.5),
     mkplusConstant(0.009),
     mHillCoefficient(5),
     mk1(0.047),
     mk2(8),
     mDt(0.01)
{
}

template<unsigned DIM>
Zuluetaforce<DIM>::~Zuluetaforce()
{
}

template<unsigned DIM>
void Zuluetaforce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    ///\todo: check whether this line influences profiling tests - if so, we should remove it.
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("Zuluetaforce is to be used with a VertexBasedCellPopulation only");
    }

    if (mRestLengths.size() == 0)
    {
        InitializeRestLengths(static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation));
        InitializeMyosinLevel(static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation));
        InitializeLineTension(static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation));
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    UpdateMyosinLevelLineTension(p_cell_population);


    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();

    /*
     * Check if a subclass of AbstractTargetAreaModifier is being used by
     * interrogating the first cell (assuming either all cells, or no cells, in
     * the population have the CellData item "target area").
     */
    bool using_target_area_modifier = p_cell_population->Begin()->GetCellData()->HasItem("target area");

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> target_areas(num_elements);
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);

        if (using_target_area_modifier)
        {
            target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
        }
        else
        {
            target_areas[elem_index] = mTargetAreaParameter;
        }
    }

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the CellPopulation, evaluated at the position of the vertex. This
         * free energy is the sum of the free energies of all CellPtrs in
         * the cell population. The free energy of each CellPtr is comprised of two
         * terms - an area deformation energy and adhesion energy term.
         *
         * Note that since the movement of this Node only affects the free energy
         * of the CellPtrs containing it, we can just consider the contributions
         * to the free energy gradient from each of these CellPtrs.
         */
        c_vector<double, DIM> area_elasticity_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> line_tension_contribution = zero_vector<double>(DIM);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned elem_index = p_element->GetIndex();
            unsigned num_nodes_elem = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            // Add the force contribution from this cell's area elasticity (note the minus sign)
            c_vector<double, DIM> element_area_gradient =
                    p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
            area_elasticity_contribution -= GetAreaElasticityParameter()*(element_areas[elem_index] -
                    target_areas[elem_index])*element_area_gradient;
            
            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);

            // Compute the line tension parameter for each of these edges - be aware that this is half of the actual
            // value for internal edges since we are looping over each of the internal edges twice

            double previous_edge_line_tension_parameter = GetLineTensionParameter(p_previous_node, p_this_node, *p_cell_population);
            double next_edge_line_tension_parameter = GetLineTensionParameter(p_this_node, p_next_node, *p_cell_population);


            // Compute the gradient of each these edges, computed at the present node
            c_vector<double, DIM> previous_edge_gradient =
                    -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
            c_vector<double, DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

            // Add the force contribution from cell-cell and cell-boundary line tension (note the minus sign)
            line_tension_contribution -= previous_edge_line_tension_parameter*previous_edge_gradient +
                    next_edge_line_tension_parameter*next_edge_gradient;
        }

        c_vector<double, DIM> force_on_node = area_elasticity_contribution + line_tension_contribution;
        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);
    }
}
template<unsigned DIM>
void Zuluetaforce<DIM>::InitializeRestLengths(VertexBasedCellPopulation<DIM>* p_cell_population)
{
    std::vector<double> rest_lengths;
    unsigned numEdges = p_cell_population->rGetMesh().GetNumEdges();
    rest_lengths.resize(numEdges);
    // Iterate over edges in the cell population
    for (unsigned edge_index=0; edge_index<numEdges; edge_index++)
    {
        // Get the rest length of this edge
        double rest_length = p_cell_population->rGetMesh().GetEdge(edge_index)->rGetLength();
        rest_lengths[edge_index] = rest_length;
    }
    mRestLengths = rest_lengths;
}

template<unsigned DIM>
void Zuluetaforce<DIM>::InitializeMyosinLevel(VertexBasedCellPopulation<DIM>* p_cell_population)
{
    std::vector<double> myosin_levels;
    unsigned numEdges = p_cell_population->rGetMesh().GetNumEdges();
    myosin_levels.resize(numEdges);
    // Iterate over edges in the cell population
    for (unsigned edge_index=0; edge_index<numEdges; edge_index++)
    {
        myosin_levels[edge_index] = 1.0;
    }
    mMyosinLevels = myosin_levels;
}

template<unsigned DIM>
void Zuluetaforce<DIM>::InitializeLineTension(VertexBasedCellPopulation<DIM>* p_cell_population)
{
    unsigned num_nodes = p_cell_population->GetNumNodes();
    mLineTension.resize(p_cell_population->rGetMesh().GetNumEdges());
    std::vector<unsigned> traversed_edges;
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned num_nodes_elem = p_element->GetNumNodes();
            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);
            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);
            std::set<unsigned> shared_elements = GetSharedElements(p_previous_node, p_this_node, *p_cell_population);
            unsigned edgeLocalIndex = GetEdgeLocalIndex(p_previous_node, p_this_node, *p_cell_population, shared_elements);
            double edgeLength = norm_2(p_previous_node->rGetLocation() - p_this_node->rGetLocation());

            if (std::find(traversed_edges.begin(), traversed_edges.end(), edgeLocalIndex) == traversed_edges.end())
            {
                traversed_edges.push_back(edgeLocalIndex);
                CalculateLineTension(edgeLocalIndex, edgeLength, shared_elements);
            }
        }
    }
}


template<unsigned DIM>
std::set<unsigned> Zuluetaforce<DIM>::GetSharedElements(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));
    assert(!shared_elements.empty());
    return shared_elements;
}

template<unsigned DIM>
unsigned Zuluetaforce<DIM>::GetEdgeLocalIndex(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation, std::set<unsigned> shared_elements)
{
    std::set<unsigned>::iterator iter = shared_elements.begin();
    // Get this element, its index and its number of nodes
    VertexElement<DIM, DIM>* p_element = rVertexCellPopulation.GetElement(*iter);
    unsigned num_edges_elem = p_element->GetNumEdges();
    unsigned edgeLocalIndex = 0;
    for (unsigned i = 0; i < num_edges_elem; i++)
    {
        if (p_element->GetEdge(i)->ContainsNode(pNodeA) && p_element->GetEdge(i)->ContainsNode(pNodeB))
        {
            edgeLocalIndex = p_element->GetEdge(i)->GetIndex();
            break;
        }
    }
    return edgeLocalIndex;
}

template<unsigned DIM>
double Zuluetaforce<DIM>::GetLineTensionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    std::set<unsigned> shared_elements = GetSharedElements(pNodeA, pNodeB, rVertexCellPopulation);
    unsigned edgeLocalIndex = GetEdgeLocalIndex(pNodeA, pNodeB, rVertexCellPopulation, shared_elements);
    double line_tension_parameter_in_calculation = mLineTension[edgeLocalIndex];
    if (shared_elements.size() == 2)
    {
        line_tension_parameter_in_calculation *= 0.5;
    }
    return line_tension_parameter_in_calculation;
}

template<unsigned DIM>
void Zuluetaforce<DIM>::UpdateMyosinLevelLineTension(VertexBasedCellPopulation<DIM>* p_cell_population)
{
    unsigned num_nodes = p_cell_population->GetNumNodes();
    std::vector<unsigned> traversed_edges;
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned num_nodes_elem = p_element->GetNumNodes();
            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);
            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);
            std::set<unsigned> shared_elements = GetSharedElements(p_previous_node, p_this_node, *p_cell_population);
            unsigned edgeLocalIndex = GetEdgeLocalIndex(p_previous_node, p_this_node, *p_cell_population, shared_elements);
            double edgeLength = norm_2(p_previous_node->rGetLocation() - p_this_node->rGetLocation());

            if (std::find(traversed_edges.begin(), traversed_edges.end(), edgeLocalIndex) == traversed_edges.end())
            {
                traversed_edges.push_back(edgeLocalIndex);
                CalculateMyosinLevel(edgeLocalIndex, edgeLength);
                CalculateLineTension(edgeLocalIndex, edgeLength, shared_elements);
            }
        }
    }
}

template<unsigned DIM>
void Zuluetaforce<DIM>::CalculateLineTension(unsigned edgeLocalIndex, double edgeLength, std::set<unsigned> shared_elements)
{
    double line_tension_in_calculation = mLineTension[edgeLocalIndex];
    if (shared_elements.size() == 2)
    {
        line_tension_in_calculation = GetNonBoundaryLineTensionParameter(edgeLocalIndex, edgeLength);
    }
    else
    {
        line_tension_in_calculation = GetBoundaryLineTensionParameter(edgeLocalIndex, edgeLength);
    }
    mLineTension[edgeLocalIndex] = line_tension_in_calculation;
}

template<unsigned DIM>
void Zuluetaforce<DIM>::CalculateMyosinLevel(unsigned edgeLocalIndex, double edgeLength)
{

    double myosin_level_in_calculation = mMyosinLevels[edgeLocalIndex];
    double deformation = edgeLength - mRestLengths[edgeLocalIndex];
    double step_function = deformation > 0 ? 1 : 0;

    double kplus = mkplusConstant + step_function * mAmax * pow(deformation, mHillCoefficient) / (pow(deformation, mHillCoefficient) + pow(mK, mHillCoefficient));
    double kminus = mk1 * exp(-mk2 * mLineTension[edgeLocalIndex]);
    myosin_level_in_calculation += mDt * (kplus * mmc - kminus * myosin_level_in_calculation);
    
    mMyosinLevels[edgeLocalIndex] = myosin_level_in_calculation;
}


template<unsigned DIM>
double Zuluetaforce<DIM>::GetBoundaryLineTensionParameter(unsigned edgeIndex, double edgeLength)
{
    double deformation = edgeLength - mRestLengths[edgeIndex];
    double deformation_contribution = mmu * deformation;
    double myosin_contribution = mbeta * mMyosinLevels[edgeIndex];
    return myosin_contribution + deformation_contribution;
}

template<unsigned DIM>
double Zuluetaforce<DIM>::GetNonBoundaryLineTensionParameter(unsigned edgeIndex, double edgeLength)
{
    return mgamma + mmu * (edgeLength - mRestLengths[edgeIndex]);
}

template<unsigned DIM>
double Zuluetaforce<DIM>::GetAreaElasticityParameter()
{
    return mAreaElasticityParameter;
}

template<unsigned DIM>
double Zuluetaforce<DIM>::GetTargetAreaParameter()
{
    return mTargetAreaParameter;
}

template<unsigned DIM>
void Zuluetaforce<DIM>::SetAreaElasticityParameter(double areaElasticityParameter)
{
    mAreaElasticityParameter = areaElasticityParameter;
}

template<unsigned DIM>
void Zuluetaforce<DIM>::SetTargetAreaParameter(double targetAreaParameter)
{
    mTargetAreaParameter = targetAreaParameter;
}

template<unsigned DIM>
void Zuluetaforce<DIM>::SetTensionParameter(double beta, double mu, double gamma)
{
    mbeta = beta;
    mmu = mu;
    mgamma = gamma;
}

template<unsigned DIM>
void Zuluetaforce<DIM>::Setmc(double mc)
{
    mmc = mc;
}

template<unsigned DIM>
void Zuluetaforce<DIM>::Setkplus(double Amax, double K, double kplusConstant, double HillCoefficient)
{
    mAmax = Amax;
    mK = K;
    mkplusConstant = kplusConstant;
    mHillCoefficient = HillCoefficient;
}

template<unsigned DIM>
void Zuluetaforce<DIM>::Setkminus(double k1, double k2)
{
    mk1 = k1;
    mk2 = k2;
}

template<unsigned DIM>
void Zuluetaforce<DIM>::SetDt(double dt)
{
    mDt = dt;
}


template<unsigned DIM>
void Zuluetaforce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AreaElasticityParameter>" << mAreaElasticityParameter << "</AreaElasticityParameter>\n";
    *rParamsFile << "\t\t\t<TargetAreaParameter>" << mTargetAreaParameter << "</TargetAreaParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class Zuluetaforce<1>;
template class Zuluetaforce<2>;
template class Zuluetaforce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Zuluetaforce)

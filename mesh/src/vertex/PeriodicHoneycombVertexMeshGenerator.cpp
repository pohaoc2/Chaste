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

#include "PeriodicHoneycombVertexMeshGenerator.hpp"
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

PeriodicHoneycombVertexMeshGenerator::PeriodicHoneycombVertexMeshGenerator(AbstractMeshReader<2, 2>& rMeshReader,
                                                           double cellRearrangementThreshold,
                                                           double t2Threshold)
{
    assert(cellRearrangementThreshold > 0.0);
    assert(t2Threshold > 0.0);

    std::vector<Node<2>*> nodes;
    std::vector<VertexElement<2,2>*> elements;
    std::vector<double> node_data;
    unsigned num_nodes = rMeshReader.GetNumNodes();
    unsigned num_elements = rMeshReader.GetNumElements();


    for (unsigned i = 0; i < num_nodes; i++)
    {
        node_data = rMeshReader.GetNextNode();
        unsigned is_boundary_node = (bool)node_data[2];
        node_data.pop_back();
        nodes.push_back(new Node<2>(i, node_data, is_boundary_node));
    }
    rMeshReader.Reset();

    for (unsigned elem_index = 0; elem_index < num_elements; elem_index++)
    {
        // Get the data for this element
        ElementData element_data = rMeshReader.GetNextElementData();

        // Get the nodes owned by this element
        std::vector<Node<2>*> element_nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j = 0; j < num_nodes_in_element; j++)
        {
            element_nodes.push_back(nodes[element_data.NodeIndices[j]]);
        }

        // Use nodes and index to construct this element
        VertexElement<2, 2>* p_element = new VertexElement<2, 2>(elem_index, element_nodes);
        elements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = (unsigned)element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }
    GenerateEdgesFromElements(elements);
    unsinged numElementsAcross = 4;
    mpMesh = boost::make_shared<Periodic2dVertexMesh>(numElementsAcross, nodes, elements, cellRearrangementThreshold, t2Threshold);
}

boost::shared_ptr<MutableVertexMesh<2,2> > PeriodicHoneycombVertexMeshGenerator::GetMesh()
{
    EXCEPTION("A cylindrical mesh was created but a normal mesh is being requested.");
    return mpMesh; // Not really
}

boost::shared_ptr<Periodic2dVertexMesh> PeriodicHoneycombVertexMeshGenerator::GetCylindricalMesh()
{
    return boost::static_pointer_cast<Periodic2dVertexMesh>(mpMesh);
}

/*
 *
 * Solution: solution data structre
 *
 *  Created on: 16/09/2015
 *      Author: bruno
 */

#ifndef SOLUTION_H_
#define SOLUTION_H_
#include <assert.h>
#include <vector>
#include <string>
#include <graph_access.h>

class Solution
{

public:

	Solution(graph_access *G);

	// add a vertex to the solution

	void addVertex(const NodeID v);

	// remove a vertex from the solution

	void removeVertex(const NodeID v);

	// randomly add a free vertex to the solution

	void addRandomVertex();

	// find a (1,2)-swap improvment

	bool twoImprovement();
	bool candTwoImprovement(NodeID candidate);

	// find a (omega,1)-swap improvment

	bool omegaImprovement();
	bool candOmegaImprovement(NodeID candidate);

	// randomly insert k vertices into the solution

	void force(NodeID k);

	void force_candidate(NodeID cand);

	// check if the solution is maximal

	bool isMaximal()
	{
		return free_size_ == 0;
	}

	// perform integrity check on the solution. returns true if
	// the solution is ok, and false, otherwise. for testing purposes
	// only.

	bool integrityCheck() const;

	// return the current solution weight

	NodeWeight weight() const
	{
		return weight_;
	}

	// return the current size of the independent set
	
	NodeID size() const
	{
		return solution_size_;
	}	

	// return a vector indicating the state of each vertex
	
	std::vector<NodeID> solution() const
	{
		std::vector<NodeID> sol;
		for (NodeID idx = 0; idx < G->number_of_nodes(); idx++) {
			if(position_[idx] < solution_size_) {
				sol.push_back(1);
			} else {
				sol.push_back(0);
			}
		}
		return sol;
	}

	// return a vector with the vertices in the solution

	std::vector<NodeID> i_set() const
	{
		std::vector<NodeID> iset;
		for (NodeID idx = 0; idx < G->number_of_nodes(); idx++) {
			if(position_[idx] < solution_size_) 
				iset.push_back(idx);
		}
		return iset;
	}


private:

	// problem instance

	graph_access *G;

	// the solution_ vector is partitioned into three blocks: first vertices in the solution, then 
	// the free vertices (i.e., vertices that are not adjacent to any vertex in the solution), and 
	// finally the non-solution vertices that are not free

	std::vector<NodeID> solution_;

	// size of the solution verticies partition

	NodeID solution_size_;

	// size of the free vertices partition

	NodeID free_size_;

	// for each vertex, the number of adjacent vertices that are on the solution

	std::vector<NodeID> tightness_;

	// position of each vertex in the solution_ vector

	std::vector<NodeID> position_;

	// weight of each vertex i minus the sum of the weights of its neighbors that
	// are in the independent set

	std::vector<long int> mu_;
	
	// current independent vertex weight
	
	NodeWeight weight_;

	// move a vertex from the free partition to solution partition

	void moveFreeToSolutionPartition(const NodeID v);

	// move a vertex from the free patition to non free partition

	void moveFreeToNonFreePartition(const NodeID v);

	// move a vertex from the solution partition to free partition

	void moveSolutionToFreePartition(const NodeID v);

	// move a vertex from the non free partition to free partition

	void moveNonFreeToFreePartition(const NodeID v);

}; // class Solution

#endif // #ifndef SOLUTION_H_

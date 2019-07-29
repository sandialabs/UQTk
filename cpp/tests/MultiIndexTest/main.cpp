/* =====================================================================================
                     The UQ Toolkit (UQTk) version @UQTKVERSION@
                     Copyright (@UQTKYEAR@) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.

     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
#include <iostream>
#include "math.h"
#include "Array1D.h"
#include "Array2D.h"
#include "probability.h"
#include "arrayio.h"
#include "arraytools.h"
#include "multiindex.h"
#include "assert.h"
#include "gen_defs.h"

using namespace std; 

/*************************************************
Begin main code
*************************************************/
int main(int argc, char ** argv){

	cout << "====> Testing multiindex generation:" << endl;
	int ndim = 3;
	int order = 4;
	Array2D<int> mindex;

	// Computes multiindex
	int npc = computeMultiIndex(ndim,order, mindex);
	// Basic sanity check
	CHECKEQ(npc,mindex.XSize());

	// Print multiindex set
	for (int ip=0; ip<mindex.XSize(); ip++){
		cout << ip << " : " ;
		for (int id=0; id<ndim; id++){
			cout << mindex(ip,id) << " ";
		}
		Array1D<int> mi;
		getRow(mindex,ip,mi);
		// Compute 'inverse' multiindex, i.e. given a multiindex, find its location
		int ip_inv=get_invmindex(mi);
		cout << " : " << ip_inv << endl;
		// Make sure inverse works correctly
		CHECKEQ(ip,ip_inv);
	}


	// Increase multiindex by an order
	cout << "Growing an order" << endl;

	Array2D<int> new_mindex;
	upOrder(mindex,new_mindex);

	// Compute orders of each term
	Array1D<int> orders;
	getOrders(new_mindex,orders);

	// Print the new multiindex
	for (int ip=0; ip<new_mindex.XSize(); ip++){
		cout << ip << " : " ;
		for (int id=0; id<ndim; id++){
			cout << new_mindex(ip,id) << " ";
		}
		cout << "; degree " << orders(ip) << endl;
		// Check that orders grew
		CHECKEQ(orders(ip),(new_mindex(ip,0)+new_mindex(ip,1)+new_mindex(ip,2)));
	}	



	return 0; 

}

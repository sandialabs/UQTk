c$$$=====================================================================================
c$$$                     The UQ Toolkit (UQTk) version @UQTKVERSION@
c$$$                     Copyright (@UQTKYEAR@) Sandia Corporation
c$$$                     http://www.sandia.gov/UQToolkit/
c$$$
c$$$     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
c$$$     with Sandia Corporation, the U.S. Government retains certain rights in this software.
c$$$
c$$$     This file is part of The UQ Toolkit (UQTk)
c$$$
c$$$     UQTk is free software: you can redistribute it and/or modify
c$$$     it under the terms of the GNU Lesser General Public License as published by
c$$$     the Free Software Foundation, either version 3 of the License, or
c$$$     (at your option) any later version.
c$$$
c$$$     UQTk is distributed in the hope that it will be useful,
c$$$     but WITHOUT ANY WARRANTY; without even the implied warranty of
c$$$     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c$$$     GNU Lesser General Public License for more details.
c$$$
c$$$     You should have received a copy of the GNU Lesser General Public License
c$$$     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
c$$$
c$$$     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
c$$$     Sandia National Laboratories, Livermore, CA, USA
c$$$=====================================================================================
      subroutine heap_ext(n,isgn,i,j,index)

       implicit none
       integer n,isgn,i,j,index

       integer l,l1,n1

       common /hpwrk/ l,l1,n1

       if (index) 90,10,80
 10    n1 = n
       l  = 1+n/2
 20    l  = l-1
 30    l1 = l
 40    i  = l1+l1
       if (i-n1) 50,60,70
 50    j  = i+1
       index = -2
       return
 60    j  = l1
       l1 = i
       index = -1
       return
 70    if (l.gt.1) goto 20
       if (n1.eq.1) goto 110
       i  = n1
       n1 = n1-1
       j  = 1
       index = 1
       return
 80    if (index-1) 30,30,40
 90    if (index.eq.-1) goto 100
       if (isgn.lt.0) i=i+1
       goto 60
 100   if (isgn.le.0) goto 70
       index = 2
       return
 110   index = 0
       return
       end




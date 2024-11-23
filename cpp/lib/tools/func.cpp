/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.5
                          Copyright (2024) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
     retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is open source software: you can redistribute it and/or modify
     it under the terms of BSD 3-Clause License

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     BSD 3 Clause License for more details.

     You should have received a copy of the BSD 3 Clause License
     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file func.cpp
/// \brief Implements several functions of form \f$y=f(\lambda;x)\f$

#include <math.h>
#include <cfloat>
#include <assert.h>

#include "func.h"

#include "gen_defs.h"
#include "PCSet.h"
#include "error_handlers.h"
#include "arrayio.h"
#include "arraytools.h"


// Proportionality function f(p,x) = p * x
Array2D<double> Func_Prop(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	assert(pdim==1);
	assert(xdim==1);

	Array2D<double> y(np,nx);

	for(int ip=0;ip<np;ip++)
		for(int ix=0;ix<nx;ix++)
			y(ip,ix)=pp(ip,0)*x(ix,0);

	return y;
}

// Quadratic proportionality f(p,x) = p1 * x + p2 * x^2
Array2D<double> Func_PropQuad(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	assert(pdim==2);
	assert(xdim==1);

	Array2D<double> y(np,nx);

	for(int ip=0;ip<np;ip++)
		for(int ix=0;ix<nx;ix++)
			y(ip,ix)=pp(ip,0)*x(ix,0)+pp(ip,1)*x(ix,0)*x(ix,0);

	return y;
}

// Exponential function f(p,x) = exp ( p1 + p2 * x )
Array2D<double> Func_Exp(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	assert(pdim==2);
	assert(xdim==1);

	Array2D<double> y(np,nx);

	for(int ip=0;ip<np;ip++)
		for(int ix=0;ix<nx;ix++)
			y(ip,ix)=exp(pp(ip,0)+pp(ip,1)*x(ix,0));

	return y;
}

// Quadratic exponential function f(p,x) = exp ( p1 + p2 * x + p3 * x^2 )
Array2D<double> Func_ExpQuad(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	assert(pdim==3);
	assert(xdim==1);

	Array2D<double> y(np,nx);

	for(int ip=0;ip<np;ip++)
		for(int ix=0;ix<nx;ix++)
			y(ip,ix)=exp(pp(ip,0)+pp(ip,1)*x(ix,0)+pp(ip,2)*x(ix,0)*x(ix,0));

	return y;
}

// Constant function f(p,x) = p
Array2D<double> Func_Const(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	assert(pdim==1);

	Array2D<double> y(np,nx);

	for(int ip=0;ip<np;ip++)
		for(int ix=0;ix<nx;ix++)
			y(ip,ix)=pp(ip,0);

	return y;
}

// Linear function f(p,x) = p1 + p2 * x
Array2D<double> Func_Linear(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	assert(pdim==2);
	assert(xdim==1);

	Array2D<double> y(np,nx);

	for(int ip=0;ip<np;ip++)
		for(int ix=0;ix<nx;ix++)
			y(ip,ix)=pp(ip,0)+pp(ip,1)*x(ix,0);

	return y;
}

// Black-box function f(p,x), expects a script bb.x which reads p.dat and x.dat and writs y.dat
Array2D<double> Func_BB(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	Array2D<double> y;

	write_datafile(pp,"p.dat"); //np,pdim
    write_datafile(x,"x.dat"); //nx,xdim
    system("./bb.x");
    read_datafileVS(y,"y.dat"); //np,nx

    CHECKEQ(y.XSize(),np);
    CHECKEQ(y.YSize(),nx);

	return y;
}

// Heat transfer example function f(p,x) = T0 + (x * dw) / (Aw * p)
Array2D<double> Func_HT1(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	assert(pdim==1);
	assert(xdim==1);

	Array2D<double> y(np,nx);

	double To=273.;
   	double Aw=0.04;
    double dw=0.1;
	for(int ip=0;ip<np;ip++)
		for(int ix=0;ix<nx;ix++)
			y(ip,ix)=x(ix,0)*dw/(Aw*pp(ip,0))+To;

	return y;
}

// Heat transfer example function f(p,x) = p2 + (x * Q) / (Aw * p1)
Array2D<double> Func_HT2(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	assert(pdim==2);
	assert(xdim==1);

	Array2D<double> y(np,nx);

	double Aw=0.04;
    //double kw=1.0;
    double Q=20.0;
	for(int ip=0;ip<np;ip++)
		for(int ix=0;ix<nx;ix++)
		    y(ip,ix)=x(ix,0)*Q/(Aw*pp(ip,0))+pp(ip,1);

	return y;
}

// Fractional power example function f(p,x) = p1 + p2 * x + p3 * x^2 + p4 * (x + 1)^3.5
Array2D<double> Func_FracPower(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	assert(pdim==4);
	assert(xdim==1);


	Array2D<double> y(np,nx);

	for(int ip=0;ip<np;ip++)
		for(int ix=0;ix<nx;ix++)
			y(ip,ix)=pp(ip,0)+pp(ip,1)*x(ix,0)+pp(ip,2)*pow(x(ix,0),2.)+pp(ip,3)*pow(x(ix,0)+1.,3.5);

	return y;
}

// Exponential example function for quick sketch f(p,x) = p2 * exp( x * p1 ) - 2
Array2D<double> Func_ExpSketch(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	assert(pdim==2);
	assert(xdim==1);

	Array2D<double> y(np,nx);

	for(int ip=0;ip<np;ip++)
		for(int ix=0;ix<nx;ix++)
			y(ip,ix)=pp(ip,1)*exp(x(ix,0)*pp(ip,0))-2.0;

	return y;
}

// Function that returns all inputs as a vector f(p,x_i) = p_i for i=1,2,..,nx
Array2D<double> Func_Inputs(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	assert(pdim==nx);
	assert(xdim==1);

	Array2D<double> y(np,nx);

	for(int ip=0;ip<np;ip++)
		for(int ix=0;ix<nx;ix++)
			y(ip,ix)=pp(ip,ix);

	return y;
}

// LU expansion f(p,x) = \sum_k p_k \Psi_k(x)
Array2D<double> Func_PCl(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	Array2D<int> mindexx;

	string pctype="LU";
    read_datafileVS(mindexx,"mindexx.dat");

    // Size checks
    CHECKEQ(mindexx.XSize(),pdim);
    CHECKEQ(mindexx.YSize(),xdim);


	PCSet surrmodel("NISPnoq",mindexx,pctype);

    Array2D<double> pcinput=x;

    Array2D<double> y(np,nx);

    for(int ip=0;ip<np;ip++){
    	Array1D<double> cf_this;
    	getRow(pp,ip,cf_this);
        Array1D<double> moutput;
        surrmodel.EvalPCAtCustPoints(moutput,pcinput,cf_this);

        for(int ix=0;ix<nx;ix++)
            y(ip,ix)=moutput(ix);
    }

	return y;
}

// LU expansion f(p,x) = \sum_k c_k \Psi(p,x)
Array2D<double> Func_PCx(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	Array2D<int> mindexpx;
	string pctype="LU";
    read_datafileVS(mindexpx,"mindexpx.dat");

    CHECKEQ(mindexpx.YSize(),pdim+xdim);

    Array1D<double> pccfpx(mindexpx.XSize());
    read_datafile_1d(pccfpx,"pccfpx.dat");

	PCSet surrmodel("NISPnoq",mindexpx,pctype);

    Array2D<double> pcinput(nx*np,pdim+xdim);
    for(int ip=0;ip<np;ip++){
	    for(int ix=0;ix<nx;ix++){
            for(int j=0;j<pdim;j++)
                pcinput(ip+ix*np,j)=pp(ip,j);
            for(int j=0;j<xdim;j++)
                pcinput(ip+ix*np,j+pdim)=x(ix,j);
        }
    }
    Array1D<double> model_samples_planar;
    surrmodel.EvalPCAtCustPoints(model_samples_planar,pcinput,pccfpx);

	Array2D<double> y(np,nx);

    for (int ip=0;ip<np;ip++)
        for (int ix=0;ix<nx;ix++)
            y(ip,ix)=model_samples_planar(ip+ix*np);

	return y;
}

// LU expansion f(p,x_i) = \sum_k c_k \Psi_k(p_i) when the multiindex is shared for all i
Array2D<double> Func_PC(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
	int* pred_mode=(int*) funcinfo;

	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	Array2D<int> mindexp;

	string pctype="LU";
    read_datafileVS(mindexp,"mindexp.dat");
    CHECKEQ(mindexp.YSize(),pdim);

    Array2D<double> pccf_all(mindexp.XSize(),nx);
    if ((*pred_mode)==0)
	    read_datafile(pccf_all,"pccf_all.dat");
    else
		read_datafile(pccf_all,"pccf_all_pred.dat");

	PCSet surrmodel("NISPnoq",mindexp,pctype);

    Array2D<double> pcinput=pp;

    Array2D<double> y(np,nx);

    for(int ix=0;ix<nx;ix++){
        Array1D<double> pccf_this;
        getCol(pccf_all,ix,pccf_this);
        Array1D<double> moutput;
        surrmodel.EvalPCAtCustPoints(moutput,pcinput,pccf_this);

        for(int ip=0;ip<np;ip++)
            y(ip,ix)=moutput(ip);
    }


	return y;
}

// LU expansion f(p,x_i) = \sum_k c_k \Psi_k(p_i) when multiindex can be different for each i
Array2D<double> Func_PCs(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo)
{
    int* pred_mode=(int*) funcinfo;

	int np=p.XSize();
	int nx=x.XSize();
	int xdim=x.YSize();

	Array2D<double> pp=augment(p,fixindnom);

	int pdim=pp.YSize();

	string pctype="LU";

	Array1D<Array2D<int> > mindices(nx);
    Array1D<Array1D<double> > pccfs(nx);

   	for (int ix=0;ix<nx;ix++){
        char filename[25];
        int nn;

        if ((*pred_mode)==0)
            nn=sprintf(filename,"mindexp.%d.dat",ix);
        else
            nn=sprintf(filename,"mindexp.%d_pred.dat",ix);
        read_datafileVS(mindices(ix), filename);
        assert(mindices(ix).YSize()==pdim);

        if ((*pred_mode)==0)
            nn=sprintf(filename,"pccfp.%d.dat",ix);
        else
            nn=sprintf(filename,"pccfp.%d_pred.dat",ix);
        pccfs(ix).Resize(mindices(ix).XSize());

        read_datafile_1d(pccfs(ix), filename);

	}



    Array2D<double> pcinput=pp;

    Array2D<double> y(np,nx);

    for(int ix=0;ix<nx;ix++){
		PCSet surrmodel("NISPnoq",mindices(ix),pctype);

        Array1D<double> moutput;
        surrmodel.EvalPCAtCustPoints(moutput,pcinput,pccfs(ix));

        for(int ip=0;ip<np;ip++)
            y(ip,ix)=moutput(ip);
    }


	return y;
}

// Augments a parameter matrix with 'fixed' columns given indices and nominal values of those
Array2D<double> augment(Array2D<double>& p, Array2D<double>& fixindnom){

   	int fdim=fixindnom.XSize();
	if (fdim>0){
	    quicksort3(fixindnom,0, fdim-1,0);
    }


    Array2D<double> pp=p;


   	for (int i=0;i<fdim;i++){
        Array1D<double> nomCol(pp.XSize(),fixindnom(i,1));
        pp.insertCol(nomCol,(int)fixindnom(i,0));
	}


	return pp;
}


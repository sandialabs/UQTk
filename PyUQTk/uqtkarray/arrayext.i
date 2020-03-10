//=====================================================================================
//
//                      The UQ Toolkit (UQTk) version 3.1.0
//                          Copyright (2020) NTESS
//                        https://www.sandia.gov/UQToolkit/
//                        https://github.com/sandialabs/UQTk
//
//     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
//     retains certain rights in this software.
//
//     This file is part of The UQ Toolkit (UQTk)
//
//     UQTk is open source software: you can redistribute it and/or modify
//     it under the terms of BSD 3-Clause License
//
//     UQTk is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     BSD 3 Clause License for more details.
//
//     You should have received a copy of the BSD 3 Clause License
//     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
//
//     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
//     Sandia National Laboratories, Livermore, CA, USA
//=====================================================================================
/*************************************************************
// Templates for Array1D and 2D types
*************************************************************/
%template(intArray1D) Array1D<int>;
%template(dblArray1D) Array1D<double>;
%template(strArray1D) Array1D<string>;
%template(intArray2D) Array2D<int>;
%template(dblArray2D) Array2D<double>;

/*************************************************************
// Extend Array1D class for easy python use
*************************************************************/

// Array 1D
%define Array1DExtend(name, T)
%extend name{
    T __getitem__(int index) {
        return (*self)[index];
    }
    Array1D<T> __getitem__(PyObject *slice) {
        Py_ssize_t start, stop, step;
        Py_ssize_t length = (*self).Length();
        Py_ssize_t slicelength;
        #ifdef PYTHON3
            PySlice_GetIndicesEx((PyObject*)slice,length,&start,&stop,&step,&slicelength);
        #endif
        #ifdef PYTHON2
            PySlice_GetIndicesEx((PySliceObject*)slice,length,&start,&stop,&step,&slicelength);
        #endif

        Array1D<T> vnew(slicelength);
        int place = 0;
        for (int i = start; i < stop; i=i+step){
            vnew(place) = (*self)[i];
            place += 1;
        }
        return vnew;
    }
    int __len__(){
        return (*self).Length();
    }
    void __setitem__(int i, T j){
        (*self)[i] = j;
    }
    void __setitem__(vector<int> index, vector<T> vin){
        // multiple index items to in vector at one time
        // In python, both index and in must be lists
        for (int i = 0; i < index.size(); i++){
            (*self)[index[i]] = vin[i];
        }
    }
    Array1D<T> __mul__(T a){
        int l = (*self).Length();
        Array1D<T> newArray(l,0);
        for (int i = 0; i < l; i++){
            newArray[i] = a*(*self)[i];
        }
        return newArray;
    }
    Array1D<T> __rmul__(T a){
        int l = (*self).Length();
        Array1D<T> newArray(l,0);
        for (int i = 0; i < l; i++){
            newArray[i] = a*(*self)[i];
        }
        return newArray;
    }
    Array1D<T> __add__(Array1D<T> y){
        int l = (*self).Length();
        Array1D<T> newArray(l,0);
        for (int i = 0; i < l; i++){
            newArray[i] = (*self)[i] + y[i];
        }
        return newArray;
    }
    Array1D<T> __add__(T y){
        int l = (*self).Length();
        Array1D<T> newArray(l,0);
        for (int i = 0; i < l; i++){
            newArray[i] = (*self)[i] + y;
        }
        return newArray;
    }
    Array1D<T> __sub__(Array1D<T> y){
        int l = (*self).Length();
        Array1D<T> newArray(l,0);
        for (int i = 0; i < l; i++){
            newArray[i] = (*self)[i] - y[i];
        }
        return newArray;
    }
    Array1D<T> __div__(Array1D<T> y){
        int l = (*self).Length();
        Array1D<T> newArray(l,0);
        for (int i = 0; i < l; i++){
            newArray[i] = (*self)[i]/y[i];
        }
        return newArray;
    }
    Array1D<T> __pow__(double p){
        int l = (*self).Length();
        Array1D<T> newArray(l,0);
        for (int i = 0; i < l; i++){
            newArray[i] = pow((*self)[i],p);
        }
        return newArray;
    }
    Array1D<T> copy(){
        int l = (*self).Length();
        Array1D<T> newArray(l,0);
        for (int i = 0; i < l; i++){
            newArray[i] = (*self)[i];
        }
        return newArray;
    }
    string __repr__(){
        stringstream ss;
        // print contents of array as strings
        int l = (*self).Length();
        ss << "Array1D<" << (*self).type() << ">(";
        ss << (*self).Length() << ")" << endl;
        int imax = 10;

        ss << "[";
        for (int i = 0; i < l-1; i++){
            // ss << setw(8) << (*self)[i] << ", ";
            ss << (*self)[i] << ", ";
            if (i == imax){
                ss << "..., ";
                break;
            }
        }
        int m = min(4,l-imax);
        if (l >= imax){
            for (int i = l-m+1; i < l-1; i++){
                ss << (*self)[i] << ", ";
            }
        }
        if (l > 0){ss << (*self)[l-1];}
        ss << "]";

        return ss.str();
    }
    vector<int> shape(){
        vector<int> s(1,0);
        s[0] = (*self).XSize();
        return s;
    }


}
%enddef

// Array 1D
%define Array1DStrExtend(name, T)
%extend name{
    T __getitem__(int index) {
        return (*self)[index];
    }
    int __len__(){
        return (*self).Length();
    }
    void __setitem__(int i, T j){
        (*self)[i] = j;
    }
    void __setitem__(vector<int> index, vector<T> vin){
        // multiple index items to in vector at one time
        // In python, both index and in must be lists
        for (int i = 0; i < index.size(); i++){
            (*self)[index[i]] = vin[i];
        }
    }
    Array1D<T> copy(){
        int l = (*self).Length();
        Array1D<T> newArray(l);
        for (int i = 0; i < l; i++){
            newArray[i] = (*self)[i];
        }
        return newArray;
    }
    string __repr__(){
        stringstream ss;
        // print contents of array as strings
        int l = (*self).Length();
        ss << "Array1D<" << (*self).type() << ">(";
        ss << (*self).Length() << ")" << endl;
        int imax = 10;

        ss << "[";
        for (int i = 0; i < l-1; i++){
            ss << (*self)[i] << ", ";
            if (i == imax){
                ss << "..., ";
                break;
            }
        }
        int m = min(4,l-imax);
        if (l >= imax){
            for (int i = l-m+1; i < l-1; i++){
                ss << (*self)[i] << ", ";
            }
        }
        if (l > 0){ss << (*self)[l-1];}
        ss << "]";

        return ss.str();
    }
    vector<int> shape(){
        vector<int> s(1,0);
        s[0] = (*self).XSize();
        return s;
    }
}
%enddef

Array1DExtend(Array1D<int>, int);
Array1DExtend(Array1D<double>, double);
Array1DStrExtend(Array1D<string>, string);

/*************************************************************
// Extend Array2D classes for easy python use
*************************************************************/

// Array2D
%define Array2DExtend(name, T)
%extend name{
    T __getitem__(vector<int> v) {
        return (*self)[v[0]][v[1]];
    }
    Array2D<T> __getitem__(PyObject *slices){
        PyObject* slice1;
        PyObject* slice2;
        slice1 = PyTuple_GetItem(slices, 0);
        slice2 = PyTuple_GetItem(slices, 1);
        //PySliceObject *s1 = (PySliceObject*)slice1; // recast pointer to proper type
        //PySliceObject *s2 = (PySliceObject*)slice2; // recast pointer to proper type

        Py_ssize_t start1 = 0, stop1 = 0, step1 = 0, slicelength1 = 0;
        Py_ssize_t start2 = 0, stop2 = 0, step2 = 0, slicelength2 = 0;
        Py_ssize_t len1 = (*self).XSize();
        Py_ssize_t len2 = (*self).YSize();

        #ifdef PYTHON3
            PySlice_GetIndicesEx(slice1,len1,&start1,&stop1,&step1,&slicelength1);
            PySlice_GetIndicesEx(slice2,len2,&start2,&stop2,&step2,&slicelength2);
        #endif

        #ifdef PYTHON2
            PySliceObject *s1 = (PySliceObject*)slice1; // recast pointer to proper type
            PySliceObject *s2 = (PySliceObject*)slice2; // recast pointer to proper type
            PySlice_GetIndicesEx(s1,len1,&start1,&stop1,&step1,&slicelength1);
            PySlice_GetIndicesEx(s2,len2,&start2,&stop2,&step2,&slicelength2);
        #endif


        Array2D<T> vnew(slicelength1,slicelength2);
        int p1 = 0, p2 = 0;
        for (int i=start1; i<stop1; i=i+step1){
            p2 = 0;
            for (int j=start2; j<stop2; j=j+step2){
                vnew(p1,p2) = (*self)[i][j];
                p2 += 1;
            }
            p1 += 1;
        }

        return vnew;
    }
    Array1D<T> __getitem__(int row) {
        (*self).getRow(row);
        return (*self).arraycopy;
    }
    int __len__(){
    	return (*self).XSize();
    }
    void __setitem__(vector<int> v, T j){
    	(*self)(v[0],v[1]) = j;
    }
    vector<int> shape(){
        vector<int> s(2,0);
        s[0] = (*self).XSize();
        s[1] = (*self).YSize();
        return s;
    }
    string __repr__(){
        stringstream ss;
        stringstream sstemp;
        // print contents of array as strings
        int lx = (*self).XSize();
        int ly = (*self).YSize();
        ss << "Array2D<" << (*self).type() << ">(";
        ss << lx << ", ";
        ss << ly << ")" << endl;

        //find # digits for number of rows, lx
        int digits = 1, pten=10;
        while ( pten <= lx ) { digits++; pten*=10; }

        // find max width (number of digits) for printing
        double test = 0.0;
        int w0 = 1, w1 = 1;
        for (int k = 0; k < lx*ly; k++){
            test = (*self).data_[k];
            sstemp.str("");
            sstemp << test;
            w1 = sstemp.str().length();
            w0 = max(w0,w1);
        }
        int w = w0;

        // size of columns for printing, dependent on width
        int imax, jmax;
        if (w >= 8){
            imax = 10;
            jmax = 8-1;
        }
        else if (w < 8){
            imax = 10;
            jmax = 12-1;
        }

        // print array
        for (int i = 0; i < lx; i++){
            // print row number
            ss << "[" << setw(digits) << i;
            ss << "]  ";
            ss << setw(2) << "[";
            for (int j=0; j<ly-1; j++){
                ss << setw(w) << (*self)[i][j] << ", ";
                if (j == jmax){
                    ss << "..." << endl;
                    ss << setw(8) << "  ";
                    break;
                }
            }
            // for each row, print after "...""
            int m = min(4,ly-jmax);
            if (ly >= jmax){
                for (int j=ly-m+1; j<ly-1; j++){
                    ss << setw(w) << (*self)[i][j] << ", ";
                }
            }
            // print last element, or first if ly == 1
            if (ly > 1){
                ss << setw(w) << (*self)[i][ly-1] << "]";
            }
            //print only if # of columns is 1
            else if (ly == 1){
                ss << setw(w) << (*self)[i][ly-1] << "]";
            }
            if (i == imax){
                ss << "\n";
                ss << setw(w) << "...," << endl;
                break;
            }
            ss << "\n";
        }

        // print last 4 rows
        int mx = min(4,lx-imax);
        if (lx >= imax){
            for (int i = lx-mx+1; i < lx; i++){
                ss << "[" << setw(digits) << i;
                ss << "]  ";
                ss << setw(2) << "[";
                for (int j=0; j<ly-1; j++){
                    ss << setw(w) << (*self)[i][j] << ", ";
                    if (j == jmax){
                        ss << "..." << endl;
                        ss << setw(8) << "  ";
                        break;
                    }
                }
                // for each row, print after ...
                int m = min(4,ly-jmax);
                if (ly >= jmax){
                    for (int j=ly-m+1; j<ly-1; j++){
                        ss << setw(w) << (*self)[i][j] << ", ";
                    }
                }
                // print last element if ly == 1
                if (ly > 1){
                    ss << setw(w) << (*self)[i][ly-1] << "]";
                }
                else if (ly == 1){
                    ss << setw(w) << (*self)[i][ly-1] << "]";
                }
                ss << "\n";
            }
        }

        return ss.str();
    }
}
%enddef

Array2DExtend(Array2D<int>, int);
Array2DExtend(Array2D<double>, double);

// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2014, Roberto P.Palomares <r.perezpalomares@gmail.com>
// All rights reserved.
#ifndef ENERGY_MODEL_H
#define ENERGY_MODEL_H

////INITIALIZATION OF AUXILIAR STUFF
void initialize_auxiliar_stuff(
      SpecificOFStuff *ofStuff, 
      OpticalFlowData *ofCore
      );

void free_auxiliar_stuff(SpecificOFStuff *ofStuff, OpticalFlowData *ofCore);

void prepare_stuff(
    SpecificOFStuff *ofStuff1,
    OpticalFlowData *ofCore1,
    SpecificOFStuff *ofStuff2,
    OpticalFlowData *ofCore2,
    float *a,
    float *b,
    int pd,
    float **out_a,
    float **out_b
    );

void of_estimation(
    SpecificOFStuff *ofStuff,
    OpticalFlowData *ofCore,
    float *ener_N,
    float *a,  //first frame
    float *b,  //second frame
    const int ii, // initial column
    const int ij, // initial row
    const int ei, // end column
    const int ej // end row
    );

void eval_functional(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore,
          float *ener_N,
          float *a,  //first frame
          float *b,  //second frame
          const int ii, // initial column
          const int ij, // initial row
          const int ei, // end column
          const int ej // end row,
          );
#endif //ENERGY_MODEL_H
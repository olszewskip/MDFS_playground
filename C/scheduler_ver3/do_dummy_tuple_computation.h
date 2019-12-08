void compute_tuple_IGs(int tuple_index[], double IGs[]) {
   for (int i = 0; i < kDim; i++)
      IGs[i] = unif(re);
}


//#ifdef TUPLE_DIMENSIONALITY == 2
//void compute_tuple_IGs(double IGs[], int col_0,
//                                     int col_1) {
//   IGs[0] = unif(re);
//   IGs[1] = unif(re);
//}
//#elif TUPLE_DIMENSIONALITY == 3
//void compute_tuple_IGs(double IGs[], int col_0,
//                                     int col_1,
//                                     int col_2) {
//   IGs[0] = unif(re);
//   IGs[1] = unif(re);
//   IGs[2] = unif(re);
//}
//#elif TUPLE_DIMENSIONALITY == 4
//void compute_tuple_IGs(double IGs[], int col_0,
//                                     int col_1,
//                                     int col_2,
//                                     int col_3) {
//   IGs[0] = unif(re);
//   IGs[1] = unif(re);
//   IGs[2] = unif(re);
//   IGs[3] = unif(re);
//}
//#elif TUPLE_DIMENSIONALITY == 5
//void compute_tuple_IGs(double IGs[], int col_0,
//                                     int col_1,
//                                     int col_2,
//                                     int col_3,
//                                     int col_4) {
//   IGs[0] = unif(re);
//   IGs[1] = unif(re);
//   IGs[2] = unif(re);
//}
//#endif


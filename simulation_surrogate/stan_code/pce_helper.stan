functions {
  matrix get_poly_order_d(vector w, int d){
    matrix[rows(w), d+1] w_p = rep_matrix(1, rows(w), d+1);
    for (i in 1:d){
      w_p[:, i+1] = w_p[:, i] .* w;
    }
    return w_p;
  }
  matrix get_PCE(matrix w_sim, int d, matrix l_poly_coeffs, array[,] int comb, int N_comb){
    int N = rows(w_sim);
    int M = cols(w_sim);
    array[M] matrix[N, d+1] poly;
    for (m in 1:M){
      matrix[N, d+1] w_sim_i_p = get_poly_order_d(w_sim[:, m], d);
      poly[m] = w_sim_i_p * l_poly_coeffs;
    }
    matrix[N, N_comb] out = rep_matrix(1., N, N_comb);
    for (i in 1:N_comb){
      for (j in 1:M){
          out[:, i] = out[:, i] .* to_vector(poly[j,:, comb[i, j]+1]);
      }
    }
    return out;
  }
}


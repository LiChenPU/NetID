#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// This script is developed to accerate NetID annotation process

// [[Rcpp::export]]
List Heterodimer_connection_core(List pgroup, double ppm) {
  
  List res;
  double ppm_tol = ppm/1E6;
  
  for(int i=0; i < pgroup.length(); i++) {
    List temp_e = as<List>(pgroup[i]);
    IntegerVector node1Vec = temp_e[0];
    IntegerVector node2Vec = temp_e[1];
    NumericVector tempe = temp_e[2];
    NumericVector tm2 = temp_e[3];
    bool condit0 = 0;
    
    double mass2 = tm2[0];
    int nnr = tempe.size();
    
    IntegerVector dfnode1Vec;
    IntegerVector dflinktypeVec;
    IntegerVector dfnode2Vec;
    std::vector<double> ppmVec;
    
    for(int r1 =0; r1 < nnr; r1++) {
      for(int r2=0; r2 < nnr; r2++) {
        double reladiff = fabs(tempe[r1] + tempe[r2] - mass2)/mass2;
        if(reladiff < ppm_tol) {
          condit0 = 1;      
          int tnode1 = node1Vec[r2];
          int linktype = node1Vec[r1];
          
          int n2 = node2Vec[0];
          
          dfnode1Vec.push_back(tnode1);
          dflinktypeVec.push_back(linktype);
          dfnode2Vec.push_back(n2);
          ppmVec.push_back(reladiff);
        }
        
      }
    }
    
    if(condit0) {
      DataFrame df = DataFrame::create(Named("node1") = dfnode1Vec,
                                       _["linktype"] = dflinktypeVec,
                                       _["node2"] = dfnode2Vec,
                                       _["mass_dif"] = ppmVec);
      res.push_back(df);
    }
  }
  
  return(res);
}

List formula_core (StringVector formulas, int sign) {
  
  List res;
  
  for(int i = 0; i < formulas.size(); i++) {
    
    String fml = formulas[i];
    fml.replace_all("D", "[2]H");
    std::string fmls = fml.get_cstring();
    
    int ende2 = fmls.length();
    std::vector<std::string> element2;
    std::vector<int> number2;
    int j = 0;
    vector<std::string> NumCharVec = {"-", ".", "0", "1",
                                      "2", "3", "4", "5", 
                                      "6", "7", "8", "9"};
    
    int b, k, m;
    
    while (j <= ende2-1) {
      
      // 1st if condition
      if(fmls.substr(j,1) == "[") {
        b = j;
        
        while (fmls.substr(j,1) != "]") {
          j++;
        }
        k = j;
        
        while (std::find(NumCharVec.begin(), NumCharVec.end(), fmls.substr(j,1)) == NumCharVec.end()) {
          j++;
        }
        m = j -1;
        
        element2.push_back(fmls.substr(b, m-b+1));
      }
      
      // 2nd if condition
      if(std::find(NumCharVec.begin(), NumCharVec.end(), fmls.substr(j,1)) == NumCharVec.end()) {
        k = j;
        while (std::find(NumCharVec.begin(), NumCharVec.end(), fmls.substr(j,1)) == NumCharVec.end()) {
          j++;
        }
        
        m = j -1;
        j--;
        element2.push_back(fmls.substr(k, m-k+1));
      }
      
      // 3rd if condition
      if(std::find(NumCharVec.begin(), NumCharVec.end(), fmls.substr(j,1)) != NumCharVec.end()) {
        k = j;
        
        while (std::find(NumCharVec.begin(), NumCharVec.end(), fmls.substr(j,1)) != NumCharVec.end()) {
          j++;
        }
        
        m = j -1;
        j--;
        
        string tmpnum = fmls.substr(k, m-k+1);
        int thisnum = stoi(tmpnum)*sign;
        
        number2.push_back(thisnum);
      }
      j++;
    }
    
    res.push_back(List::create(_["element2"] = element2, 
                               _["number2"]  = number2));
    
  }
  
  return(res);
}

bool is_negative(int element){
  return (element < 0);
}

IntegerVector eleCount_sort(IntegerVector x, StringVector y) {
  IntegerVector idx = seq_along(x) - 1;
  std::sort(idx.begin(), idx.end(), [&](int i, int j){
    String ss = y[i];
    string sss = ss.get_cstring();
    
    String xx = y[j];
    string xxx = xx.get_cstring();
    bool res = 0;
    if(sss.substr(0,1) == "[" && xxx.substr(0,1) != "[") {
      res = 1;
      return res;
    } else if (sss.substr(0,1) != "[" && xxx.substr(0,1) == "[") {
      res = 0;
      return res;
    }
    return y[i] < y[j];
  }
  );
  return x[idx];
}

StringVector element_sort (StringVector y) {
  IntegerVector idx = seq_len(y.length()) - 1;
  std::sort(idx.begin(), idx.end(), [&](int i, int j){
    String ss = y[i];
    string sss = ss.get_cstring();
    
    String xx = y[j];
    string xxx = xx.get_cstring();
    bool res = 0;
    if(sss.substr(0,1) == "[" && xxx.substr(0,1) != "[") {
      res = 1;
      return res;
    } else if (sss.substr(0,1) != "[" && xxx.substr(0,1) == "[") {
      res = 0;
      return res;
    }
    return y[i] < y[j];
  }
  );
  return y[idx];
}

// [[Rcpp::export]]
StringVector fast_calculate_formula (StringVector formular1, StringVector transformulas, int sign) {
  // This function is designed to speed up "my_calculate_formula" from lc8
  
  int signNum = sign;
  List formulaList1 = formula_core(formular1, sign = 1);
  List formulaList2 = formula_core(transformulas, sign = signNum);
  
  int nrow = formulaList1.length();
  int ncol = formulaList2.length();
  
  StringMatrix formula_mat (nrow, ncol);
  LogicalMatrix valid_mat (nrow, ncol);
  
  StringVector res;
  
  for (int r = 0; r < nrow; r++) {
    List lst1 = as<List>(formulaList1[r]);
    StringVector elements1 = lst1[0];
    IntegerVector numbers1 = lst1[1];
    
    for (int c = 0; c < ncol; c++) {
      List lst2 = as<List>(formulaList2[c]);
      StringVector elements2 = lst2[0];
      IntegerVector numbers2 = lst2[1];
      
      IntegerVector count_all = Rcpp::clone(numbers1);;
      StringVector elem_all = Rcpp::clone(elements1);
      
      int endEleIdx = elem_all.end().index();
      
      for(int e =0; e < elements2.size(); e++) {
        int matchedElement = std::find(elem_all.begin(), elem_all.end(), elements2[e]).index();
        
        if(matchedElement < endEleIdx) {
          count_all[matchedElement] = count_all[matchedElement] + numbers2[e];
        } else {
          count_all.push_back(numbers2[e]);
          elem_all.push_back(elements2[e]);
        }
      }
      
      if(std::any_of(count_all.begin(), count_all.end(), is_negative)) {
        valid_mat(r,c) = 0;
      }
      
      count_all = eleCount_sort(count_all,elem_all);
      elem_all = element_sort(elem_all);
      
      String realformula;
      
      for(int cc = 0; cc < count_all.size(); cc++){
        if(count_all[cc] == 0) {
          count_all.erase(cc);
          elem_all.erase(cc);
          cc--;
          continue;
        }
        
        int count_s = count_all[cc];
        String ssss = elem_all[cc];
        realformula = realformula.push_back(ssss).push_back(count_s);
      }
      formula_mat(r,c) = realformula;
    }
    
    //if(formular1.size() == 1 && transformulas.size() == 1) {
    //  res = formula_mat(1,1);
    //} else {
    //res = StringVector(formula_mat);
    //}
  }
  
  // for(int rr = 0; rr < nrow; rr++) {
  //   Rcout << "running here subres  <-- rr " << rr << "\n";
  //   for(int cl = 0; cl < ncol; cl++) {
  //     Rcout << "running here subres  <-- cl " << rr << " | "<< cl << "\n";
  //     String subres = formula_mat(rr,cl);
  //     res.push_back(subres);
  //   }
  // }
  res = StringVector(formula_mat);
  return(res);
}

// [[Rcpp::export]]
List propagate_heterodimer_core(DataFrame df_heterodimer, 
                                List sf,
                                StringVector propagation_category,
                                NumericVector node_mass,
                                double ppm_threshold) {
  
  List res;
  
  NumericVector parantIDVec =  df_heterodimer[5];
  NumericVector trasnformVec =  df_heterodimer[7];
  NumericVector node2Vec =  df_heterodimer[14];
  NumericVector rdbeVec = df_heterodimer[3];
  NumericVector massVec =  df_heterodimer[2];
  
  StringVector formulaVec =  df_heterodimer[1];
  StringVector parentForVec =  df_heterodimer[6];
  
  
  for(int i = 0; i < trasnformVec.size(); i++) {
    // Replace the dyplr part here
    DataFrame df = as<DataFrame>(sf[trasnformVec[i]-1]);
    IntegerVector ndidVec = df[0];
    StringVector formuVec = df[1];
    NumericVector msVec = df[2];
    NumericVector reVec = df[3];
    StringVector catVec = df[4];
    NumericVector stpsVec = df[9];
    
    StringVector us = unique(formuVec);
    IntegerVector order_row = IntegerVector::create(0);
    
    IntegerVector parentid_res, transform_res, node2_res;
    NumericVector rdbe_res, mass_res;
    StringVector parentf_res, formula_res;
    
    for(int j = 0; j < us.size(); j++) {
      for(int r = 0; r < df.nrow(); r ++) {
        if(formuVec[r] == us[j]) {
          // filter according to steps
          bool cond1 = (stpsVec[r] <= 0.01 || (stpsVec[r] >= 1 && stpsVec[r] <= 1.01));
          // filter according to category
          bool cond2 = 0;
          for(int k = 0; k < propagation_category.size(); k++) {
            if(propagation_category[k] == catVec[r]){
              cond2 = 1;
              break;
            }
          }
          // filter mass
          bool cond3 = fabs(node_mass[ndidVec[r] -1] - msVec[r]) < ppm_threshold*msVec[r];
          
          if(cond1 && cond2 && cond3) {
            // it r equals transform
            // it i equals temp
            
            order_row.push_back(r);
            
            parentid_res.push_back(parantIDVec[i]);
            transform_res.push_back(trasnformVec[i]);
            node2_res.push_back(node2Vec[i]);
            parentf_res.push_back(parentForVec[i]);
            
            mass_res.push_back(massVec[i] + msVec[r]);
            rdbe_res.push_back(rdbeVec[i] + reVec[r]);
            
            formula_res.push_back(formuVec[r]);
            break;
          }
        }
      }
    }
    
    if(order_row.size() == 1) {
      continue;
    }
    
    // Create list
    StringVector tmpStrVec = StringVector::create("");
    tmpStrVec[0] = formulaVec[i];
    formula_res = fast_calculate_formula(tmpStrVec, formula_res, 1);
    
    res.push_back(List::create(_["parent_id"]  = parentid_res, 
                               _["transform"]  = transform_res, 
                               _["node2"] = node2_res,
                               _["parent_formula"] = parentf_res,
                               _["formula"] = formula_res,
                               _["rdbe"] = rdbe_res,
                               _["mass"] = mass_res));
    
  }
  return(res);
}


int find_int_pos (int num, IntegerVector numVec) {
  for (int i = 0; i < numVec.length(); i++) {
    if(numVec[i] == num) {
      return(i);
    }
  }
  return(-1);
}

IntegerVector find_int_poss (int num, IntegerVector numVec) {
  IntegerVector res;
  for (int i = 0; i < numVec.length(); i++) {
    if(numVec[i] == num) {
      res.push_back(i);
    }
  }
  return(res);
}


void coreFUN (int i, 
              StringVector &res,
              IntegerVector &solRes, 
              IntegerVector &ilpNodeIdVec,
              CharacterVector &colNMs,
              CharacterVector &rowNMs,
              IntegerVector &core_ilpNIdVec,
              StringVector &core_annoVec,
              NumericMatrix &dis_mat,
              List &g_annotation,
              Function &f,
              IntegerVector &IP_edge_From,
              IntegerVector &IP_edge_To,
              StringVector &classVec,
              StringVector &IP_edge_cat,
              StringVector &IP_edge_linktp,
              StringVector &IP_edge_fm1,
              StringVector &IP_edge_fm2,
              IntegerVector &IP_edge_dir,
              string igraphmode) {
  
  int IlpNodeID = ilpNodeIdVec[i]; // this equals query_ilp_id

  if(std::find(colNMs.begin(), colNMs.end(), IlpNodeID) == colNMs.end()) {
    res[i] = "not exited in dm";
    return;
  }
  //std::to_string()
  if(std::find(rowNMs.begin(), rowNMs.end(), IlpNodeID) != rowNMs.end()){
    int core_pos = find_int_pos(IlpNodeID, core_ilpNIdVec);
    res[i] = core_annoVec[core_pos];
    return;
  }
  
  int tmp_pos = std::find(colNMs.begin(), colNMs.end(), IlpNodeID).index();
  NumericVector v = dis_mat(_, tmp_pos);
  double dm_min = min(v);
  double inf = std::numeric_limits<double>::infinity();
  
  if(dm_min == inf) {
    res[i] = "No edge connections.";
    return;
  }
  int np = which_min(v);
  String parent_selected = rowNMs[np];
  
  List path_nodes = f(g_annotation, 
                      Named("from", parent_selected),
                      Named("to", std::to_string(IlpNodeID)),
                      Named("mode", igraphmode),
                      Named("output", "vpath"));
  List tmp1 = as<List>(path_nodes[0]);
  IntegerVector tmp2 = as<IntegerVector>(tmp1[0]);
  StringVector tmp3 = tmp2.names();
  IntegerVector ilp_node_path(tmp3.size());
  for(int j = 0; j < tmp3.size(); j++){
    String tmp4 = tmp3[j];
    ilp_node_path[j] = stoi(tmp4.get_cstring());
  }
  
  int corePos = find_int_pos(ilp_node_path[0], core_ilpNIdVec);
  
  String temp_core = core_annoVec[corePos];
  String transform_path;
  int dir = 0;
  
  for(int k =0; k < ilp_node_path.size()-1; k++) {
    IntegerVector resFrom = find_int_poss(ilp_node_path[k], IP_edge_From);
    IntegerVector resTo = find_int_poss(ilp_node_path[k+1], IP_edge_To);
    IntegerVector resposs = intersect(resFrom, resTo);
    resposs = resposs.sort();
    
    if(resposs.size() > 0) {
      dir = IP_edge_dir[resposs[0]];
    } else if(resposs.size() == 0) {
      resFrom = find_int_poss(ilp_node_path[k], IP_edge_To);
      resTo = find_int_poss(ilp_node_path[k+1], IP_edge_From);
      resposs = intersect(resFrom, resTo);
      dir = -1;
    }
    
    string temp_sign;
    if(classVec[i] == "Artifact") {
      if(IP_edge_cat[resposs[0]] == "Oligomer") {
        temp_sign = "*";
      } else if (IP_edge_cat[resposs[0]] == "Multicharge") {
        temp_sign = "/";
      } else if (IP_edge_cat[resposs[0]] == "Heterodimer") {
        temp_sign = "+ Peak";
      } else if (dir == 1) {
        temp_sign = "+";
      } else if (dir == -1){
        temp_sign = "-";
      }
    } else {
      temp_sign = (dir == 1) ? "+" : "-";
    }
    
    String temp_linktype = IP_edge_linktp[resposs[0]];
    String temp_formula = (dir == 1) ? IP_edge_fm2[resposs[0]] : IP_edge_fm1[resposs[0]];
    
    string seperator = " ";
    transform_path =
      transform_path.get_cstring() + seperator +
      temp_sign + seperator +
      temp_linktype.get_cstring() + seperator + "-> " + 
      temp_formula.get_cstring();
  }
  string reso = transform_path.get_cstring();
  string resx= temp_core.get_cstring() + reso;
  res[i] = resx;
  
}

// [[Rcpp::export]]
StringVector path_annotate(DataFrame ilp_nodes, 
                           double solution_index,
                           double class_index,
                                DataFrame canu_met, //core_annotation_unique_*
                                DataFrame ilp_edges_anno_met, //ilp_edges_annotate_*
                                NumericMatrix dis_mat_met,
                                List g_annotation,
                                DataFrame canu_nonmet, //core_annotation_unique_*
                                DataFrame ilp_edges_anno_nonmet, //ilp_edges_annotate_*
                                NumericMatrix dis_mat_nonmet,
                                List g_anno_non) {
  
  IntegerVector solRes = as<IntegerVector>(ilp_nodes[solution_index]);
  IntegerVector ilpNodeIdVec = as<IntegerVector>(ilp_nodes[0]);
  StringVector classVec = as<StringVector>(ilp_nodes[class_index]);
  StringVector res(solRes.length());
  
  // for met
  IntegerVector core_ilpNIdVec1 = as<IntegerVector>(canu_met[0]);
  StringVector core_annoVec1 = as<StringVector>(canu_met[2]);
  
  CharacterVector colNMs1 = colnames(dis_mat_met);
  CharacterVector rowNMs1 = rownames(dis_mat_met);
  
  IntegerVector IP_edge_From1 = as<IntegerVector>(ilp_edges_anno_met[0]);
  IntegerVector IP_edge_To1 = as<IntegerVector>(ilp_edges_anno_met[1]);
  IntegerVector IP_edge_dir1 = as<IntegerVector>(ilp_edges_anno_met[8]);
  StringVector IP_edge_linktp1 = as<StringVector>(ilp_edges_anno_met[7]);
  StringVector IP_edge_fm1_1 = as<StringVector>(ilp_edges_anno_met[10]);
  StringVector IP_edge_fm2_1 = as<StringVector>(ilp_edges_anno_met[11]);
  StringVector IP_edge_cat1 = as<StringVector>(ilp_edges_anno_met[6]);
  
  // for nonmet
  IntegerVector core_ilpNIdVec2 = as<IntegerVector>(canu_nonmet[0]);
  StringVector core_annoVec2 = as<StringVector>(canu_nonmet[2]);
  
  CharacterVector colNMs2 = colnames(dis_mat_nonmet);
  CharacterVector rowNMs2 = rownames(dis_mat_nonmet);
  
  IntegerVector IP_edge_From2 = as<IntegerVector>(ilp_edges_anno_nonmet[0]);
  IntegerVector IP_edge_To2 = as<IntegerVector>(ilp_edges_anno_nonmet[1]);
  IntegerVector IP_edge_dir2 = as<IntegerVector>(ilp_edges_anno_nonmet[8]);
  StringVector IP_edge_linktp2 = as<StringVector>(ilp_edges_anno_nonmet[7]);
  StringVector IP_edge_fm1_2 = as<StringVector>(ilp_edges_anno_nonmet[10]);
  StringVector IP_edge_fm2_2 = as<StringVector>(ilp_edges_anno_nonmet[11]);
  StringVector IP_edge_cat2 = as<StringVector>(ilp_edges_anno_nonmet[6]);
  
  // igraph FUN
  Environment pkg = Environment::namespace_env("igraph");
  Function f = pkg["shortest_paths"];
  
  for(int i = 0; i < solRes.size(); i++) {
    if(solRes[i] == 0) {continue;}
    if(classVec[i] == "Metabolite" || classVec[i] == "Putative metabolite") {
      coreFUN(i, 
              res,
              solRes, 
              ilpNodeIdVec,
              colNMs1,
              rowNMs1,
              core_ilpNIdVec1,
              core_annoVec1,
              dis_mat_met,
              g_annotation,
              f,
              IP_edge_From1,
              IP_edge_To1,
              classVec,
              IP_edge_cat1,
              IP_edge_linktp1,
              IP_edge_fm1_1,
              IP_edge_fm2_1,
              IP_edge_dir1,
              "all");
    } else if (classVec[i] == "Artifact") {
      coreFUN(i, 
              res,
              solRes, 
              ilpNodeIdVec,
              colNMs2,
              rowNMs2,
              core_ilpNIdVec2,
              core_annoVec2,
              dis_mat_nonmet,
              g_anno_non,
              f,
              IP_edge_From2,
              IP_edge_To2,
              classVec,
              IP_edge_cat2,
              IP_edge_linktp2,
              IP_edge_fm1_2,
              IP_edge_fm2_2,
              IP_edge_dir2,
              "out");
    } else {
      res[i] = "Unknown";
    }
  }
  
  return res;
}

// [[Rcpp::export]]
StringVector path_annotate_met_only(DataFrame ilp_nodes, 
                           double solution_index,
                           double class_index,
                                DataFrame canu_met, //core_annotation_unique_*
                                DataFrame ilp_edges_anno_met, //ilp_edges_annotate_*
                                NumericMatrix dis_mat_met,
                                List g_annotation) {
  
  IntegerVector solRes = as<IntegerVector>(ilp_nodes[solution_index]);
  IntegerVector ilpNodeIdVec = as<IntegerVector>(ilp_nodes[0]);
  StringVector classVec = as<StringVector>(ilp_nodes[class_index]);
  StringVector res(solRes.length());
  
  // for met
  IntegerVector core_ilpNIdVec1 = as<IntegerVector>(canu_met[0]);
  StringVector core_annoVec1 = as<StringVector>(canu_met[2]);
  
  CharacterVector colNMs1 = colnames(dis_mat_met);
  CharacterVector rowNMs1 = rownames(dis_mat_met);
  
  IntegerVector IP_edge_From1 = as<IntegerVector>(ilp_edges_anno_met[0]);
  IntegerVector IP_edge_To1 = as<IntegerVector>(ilp_edges_anno_met[1]);
  IntegerVector IP_edge_dir1 = as<IntegerVector>(ilp_edges_anno_met[8]);
  StringVector IP_edge_linktp1 = as<StringVector>(ilp_edges_anno_met[7]);
  StringVector IP_edge_fm1_1 = as<StringVector>(ilp_edges_anno_met[10]);
  StringVector IP_edge_fm2_1 = as<StringVector>(ilp_edges_anno_met[11]);
  StringVector IP_edge_cat1 = as<StringVector>(ilp_edges_anno_met[6]);
  
    
  // igraph FUN
  Environment pkg = Environment::namespace_env("igraph");
  Function f = pkg["shortest_paths"];
  
  for(int i = 0; i < solRes.size(); i++) {
    if(solRes[i] == 0) {continue;}
    if(classVec[i] == "Metabolite" || classVec[i] == "Putative metabolite") {
      coreFUN(i, 
              res,
              solRes, 
              ilpNodeIdVec,
              colNMs1,
              rowNMs1,
              core_ilpNIdVec1,
              core_annoVec1,
              dis_mat_met,
              g_annotation,
              f,
              IP_edge_From1,
              IP_edge_To1,
              classVec,
              IP_edge_cat1,
              IP_edge_linktp1,
              IP_edge_fm1_1,
              IP_edge_fm2_1,
              IP_edge_dir1,
              "all");
    } else if (classVec[i] == "Artifact") {
      continue;
    } else {
      res[i] = "Unknown";
    }
  }
  
  return res;
}



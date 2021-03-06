findpolicy_Milad2 = function(initial, file = 'output.policy'){
  library(xml2)
  library(purrr)
  ## Extract vectors from XML
  xml <- xml2::read_xml(file)
  vectors <- xml2::xml_find_all(xml, "//Vector")
  get_vector <- function(v) as.numeric(strsplit(as.character(xml2::xml_contents(v)), " ")[[1]]) 
  alpha <- vectors %>% purrr::map(get_vector)
  alpha_action <- vectors %>% purrr::map_dbl(function(v) as.numeric(xml2::xml_attr(v, "action")))
  
  ## Compute dot product with initial
  a <- alpha %>% map_dbl(function(x) initial %*% matrix(x, ncol=1))

  ## Return policy of the vector which has the biggest inner product
  #alpha_action[which.max(a)]
  
  output = list(max(a),alpha_action[which.max(a)],alpha,alpha_action)
}

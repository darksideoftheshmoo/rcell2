`%|%` <- function(l1, l2){
  l1[names(l2)] <- l2
  return(l1)
}

`%|=%` <- function(l1, l2){
  common.names <- names(l2)[names(l2) %in% names(l1)]
  l1[common.names] <- l2[common.names]
  return(l1)
}

lista1 <- list(hola = "HOLA", chau = "BYE", hmm = "HMM")

lista2 <- list(hola = "HOLA", chau = "CHAUCHAS", quetal = "COMO VA")

lista1
lista2
lista1 %|% lista2

lista1 %|=% lista2

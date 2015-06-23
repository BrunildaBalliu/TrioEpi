TCfunction=function(father, mother, child){
  n=nrow(child);
  M=ncol(child);
  
  # Observed genotype frequencies
  obs = matrix(0,nrow=3,ncol=3);  
  for (i in 0:2) {
    for (j in 0:2) {
      obs[i+1,j+1] = sum(child[,1]==i&child[,2]==j)
    }
  };
  
  # Expected genotype frequencies
  p_0 = matrix(0,nrow=n,ncol=M); # = P(g^o = 0 | g^f, g^m)
  p_1 = matrix(0,nrow=n,ncol=M); # = P(g^o = 1 | g^f, g^m)
  p_2 = matrix(0,nrow=n,ncol=M); # = P(g^o = 2 | g^f, g^m)
  
  for (i in 1:n) {
    p_0[i,] = (1*(father[i,]==0)+0.5*(father[i,]==1)+0*(father[i,]==2))*(1*(mother[i,]==0)+0.5*(mother[i,]==1)+0*(mother[i,]==2))    
    p_1[i,] = 0*(father[i,]==0&mother[i,]==0)+0.5*(father[i,]==0&mother[i,]==1)+1*(father[i,]==0&mother[i,]==2)+0.5*(father[i,]==1&mother[i,]==0)+0.5*(father[i,]==1&mother[i,]==1)+0.5*(father[i,]==1&mother[i,]==2)+ 1*(father[i,]==2&mother[i,]==0)+0.5*(father[i,]==2&mother[i,]==1)+0*(father[i,]==2&mother[i,]==2)    
    p_2[i,] = (0*(father[i,]==0)+0.5*(father[i,]==1)+1*(father[i,]==2))*(0*(mother[i,]==0)+0.5*(mother[i,]==1)+1*(mother[i,]==2))
  };
  
  p = list();
  p[[1]] = p_0;
  p[[2]] = p_1;
  p[[3]] = p_2;
  
  exp = matrix(0,nrow=3,ncol=3);
  for (i in 0:2) {
    for (j in 0:2) {
      exp[i+1,j+1] = sum(p[[i+1]][,1]*p[[j+1]][,2])
    }
  };
  
  # TC test statistic and p-value
  rowc = apply(exp,2,sum);
  mrow = rowc%*%c(0,1,2)/sum(rowc);
  
  colc = apply(exp,1,sum);
  mcol = colc%*%c(0,1,2)/sum(colc);
  
  sd1 = sd2 = cor = 0;
  
  for (i in 1:3) {
    for (j in 1:3) {
      cor = cor+obs[i,j]*((i-1)-mcol)*((j-1)-mrow)
      sd1 = sd1 + (obs[i,j]*((i-1)-mcol)^2)
      sd2 = sd2 + (obs[i,j]*((j-1)-mrow)^2)
    }
  };
  
  cor = cor/(sqrt(sd1)*sqrt(sd2));
  TestStatistic = cor^2*n;
  Pvalue=pchisq(q = TestStatistic, df = 1, lower.tail = F);  
  return(list(TestStatistic=TestStatistic, Pvalue=Pvalue))
}
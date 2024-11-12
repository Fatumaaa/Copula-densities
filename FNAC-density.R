library(copula)

d <- 4
tau.values <- seq(from = 1, to = 0, length = d + 1)[-c(1,(d + 1))]
theta.values <- copGumbel@iTau(tau.values)

generate.data <- function(nr.times, n, copula.1){
  FNAC_Gumbel <- onacopula(copula.1, C(theta.values[3], 4, C(theta.values[2], 3, C(theta.values[1], c(2, 1)))))
  return(rnacopula(n = n, FNAC_Gumbel))
}

# Get one dataset
data.Gumbel <- lapply(1, generate.data, n = 50, copula.1 = "Gumbel")

theta.function <- function(d){
  ret <- matrix(ncol = (d - 1), nrow = (d - 1))
  for(i in 1:(d-1)){
    ret[i,i:(d-1)] <- i:(d-1)
  }
  return(ret)
}

matrix.theta <- theta.function(d)

n.deriv.Gumbel <- function(n, t, theta){
  s1 <- c()
  for(r in 1:n){
    for(s in 1:r){
      vec <- (s/theta)-(0:(n-1))
      s1 <- c(s1, ((-1)^s/factorial(r))*choose(r, s)*t^(r/theta)*prod(vec))
    }
  }
  return(exp(-t^(1/theta))*t^(-n)*sum(s1))
}

psi.ij.q <- function(n, q, t, theta.1, theta.2){
  vec <- (q*theta.1/theta.2) - (0:(n-1))
  return(t^(-n + q*theta.1/theta.2)*prod(vec))
}

constante.function <- function(p, q){
  return(((-1)^(p + q)/factorial(p))*choose(p, q))
}

psi.Gumbel <- function(t, theta){
  return(exp(-t^(1/theta)))
}

psi.inv.Gumbel <- function(t, theta){
  return((-log(t))^(theta))
}

deriv.inv.Gumbel <- function(t, theta){
  return(theta*(-log(t))^theta/(t*log(t)))
}

expand.grid.function <- function(i){
  eg <- expand.grid(1:i, 1:i)
  return(eg[which(eg$Var2 <= eg$Var1),])
}

list.elem <- lapply(1:(d-1), expand.grid.function)

eg <- expand.grid(as.list(data.frame(t(list.elem[[1]]))), as.list(data.frame(t(list.elem[[2]]))))
unlist.function <- function(i, eg){
  return(unlist(eg[i,]))
}
eg <- lapply(1:nrow(eg), unlist.function, eg = eg)
for(i in 3:(d-1)){
  eg.1 <- expand.grid(eg, as.list(data.frame(t(list.elem[[i]]))))
  eg <- lapply(1:nrow(eg.1), unlist.function, eg = eg.1)
}

p.q.matrix <- t(sapply(1:nrow(eg.1), unlist.function, eg = eg.1))
p.matrix <- data.frame(p.q.matrix[,seq(1,by=2, len=d-1)])
q.matrix <- data.frame(p.q.matrix[,seq(2,by=2, len=d-1)])

loglikelihood.function <- function(par,U, d){
  theta <- c()
  for(i in 1:nrow(matrix.theta)){
    theta <- c(theta, sum(exp(par[matrix.theta[i,which(!is.na(matrix.theta[i,]))]])))
  }
  
  theta.theta <- c(theta[1], theta)

  ll.individual <- function(k){
    deriv.prod.G <- function(i){
      return(deriv.inv.Gumbel(U[k,i], theta.theta[i]))    
    }
    pp <- prod(deriv.prod.G(1:d))
    
    ## From the inside out create the copula
    Cop.Dim <- function(d){
      if(d == 1){
        return(U[k,1])
      }else{
        Cop.Gumbel <- U[k,1]
        for(i in 1:(d-1)){
          Cop.Gumbel <- psi.Gumbel(psi.inv.Gumbel(U[k,i + 1], theta[i]) + psi.inv.Gumbel(Cop.Gumbel, theta[i]), theta[i]) 
        }
        return(Cop.Gumbel)
      }
    }
    
    deriv.psi <- function(p){
      return(n.deriv.Gumbel(p + 1,psi.inv.Gumbel(U[k,d], theta = theta[d-1]) + psi.inv.Gumbel(Cop.Dim(d-1), theta = theta[d-1]), theta = theta[d-1]))
    }
    
    deriv.psi.last <- sapply(p.matrix[,d-1], deriv.psi)
    
    middle.function <- function(i){
      middle <- c()
      for(j in 1:(d-1)){
        if(j == 1){
          middle <- c(middle, constante.function(p.matrix[i,j], q.matrix[i,j])*psi.inv.Gumbel(Cop.Dim(j), theta[j])^(p.matrix[i,j] - q.matrix[i,j]))
        }else{
          middle <- c(middle, constante.function(p.matrix[i,j], q.matrix[i,j])*psi.inv.Gumbel(Cop.Dim(j), theta[j])^(p.matrix[i,j] - q.matrix[i,j])*psi.ij.q(p.matrix[i,j-1] + 1, q.matrix[i,j], psi.inv.Gumbel(U[k,j], theta[j - 1]) + psi.inv.Gumbel(Cop.Dim(j-1), theta[j - 1]), theta.1 = theta[j], theta.2 = theta[j - 1]))
        }}
      
      return(prod(middle))
    }
    middle.part <- sapply(1:nrow(p.matrix), middle.function)
    ll.individual <- pp*sum(middle.part*deriv.psi.last)
    return(log(ll.individual))
  }
  inside.sum <- sapply(1:nrow(U), ll.individual)
  return(-sum(inside.sum))
}

optimize.function <- function(i){
  res.1 <- copGumbel@iTau(cor.test(data.Gumbel[[i]][,1],data.Gumbel[[i]][,2], method="kendall")[4]$estimate)
  res.2 <- copGumbel@iTau(cor.test(data.Gumbel[[i]][,1],data.Gumbel[[i]][,3], method="kendall")[4]$estimate)
  res.3 <- copGumbel@iTau(cor.test(data.Gumbel[[i]][,1],data.Gumbel[[i]][,4], method="kendall")[4]$estimate)
  initial.value <- c(log(res.1 - res.2),log(res.2 - res.3),log(res.3))
  initial.value[which(!is.finite(initial.value))] <- runif(length(which(!is.finite(initial.value))))
  op <- try(optim(par = initial.value, loglikelihood.function, U = data.Gumbel[[i]], d = d, method = "BFGS"), silent = FALSE)
  if(inherits(op, 'try-error')){
    print("Error")
    op <- NA
  }
  return(op)
}

optimize.all.data <- lapply(1:length(data.Gumbel), optimize.function)

# Get back the parameters, names it tau but that is not correct. 
theta.1 <- exp(optimize.all.data[[1]]$par[1]) + exp(optimize.all.data[[1]]$par[2]) + exp(optimize.all.data[[1]]$par[3])
theta.2 <- exp(optimize.all.data[[1]]$par[2]) + exp(optimize.all.data[[1]]$par[3]) 
theta.3 <- exp(optimize.all.data[[1]]$par[3])
print(c(theta.1, theta.2, theta.3))

## true theta values
theta.values

## Expressions for the Clayton copula
#
# n.deriv.Clayton <- function(n, t, theta){
#   s1 <- c()
#   for(r in 1:n){
#     for(s in 1:r){
#       vec <- (s/theta)-(0:(n-1))
#       #print(vec)
#       s1 <- c(s1, ((-1)^s/factorial(r))*choose(r, s)*t^(r/theta)*prod(vec))
#       #print(sum(s1))
#     }
#   }
#   return(exp(-t^(1/theta))*t^(-n)*sum(s1))
# }
# 
# n.deriv.Clayton <- function(n, t, theta){
#   a <- (-1)^(n)*(1 + theta*t)^(-1/theta - n)
#   vec <- (0:(n-1)) + 1/theta
#   return(a*prod(vec)*theta^n)
# }
# 
# psi.ij.q <- function(n, q, t, theta.1, theta.2){
#   vec <- c()
#   for(r in 0:q){
#     vec.1 <- r*theta.1/theta.2 - (0:(n-1))
#     vec <- c(vec, choose(q,r)*(-1)^(q - r)*theta.2^n*theta.1^(-q)*(1 + theta.2*t)^(r*theta.1/theta.2 - n)*prod(vec.1))
#   }
#   return(sum(vec))
# }
# 
# psi.Clayton <- function(t, theta){
#   return((1 + theta*t)^(-1/theta))
# }
# 
# psi.inv.Clayton <- function(t, theta){
#   return((1/theta)*(t^(-theta) - 1))
# }
# 
# deriv.inv.Clayton <- function(t, theta){
#   return(-t^(-theta - 1))
# }




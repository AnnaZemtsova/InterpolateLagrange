
setClass(
  "Interpolation",
  representation (y="vector",x="vector",l="list",tmpList="list"),
  prototype(y=0,x=0)
)

setGeneric(
  "symbolQuotient",
  function(object,j,i){
    standardGeneric("symbolQuotient")
  }
)

setGeneric(
  "multiply",
  function(object){
    standardGeneric("multiply")
  }
)

setGeneric(
  "getLagranjStringPolynom",
  function(object,a){
    setGeneric("getLagranjStringPolynom")
  }
)

setGeneric(
  "findLagranjPolynom",
  function(object){
    standardGeneric("findLagranjPolynom")
  }
)

setGeneric(
  "getLagranjPolynom",
  function(object){
    standardGeneric("getLagranjPolynom")
  }
)

setGeneric(
  "findBasicPolynom",
  function(object,j){
    standardGeneric("findBasicPolynom")
  }
)

setGeneric(
  "reduction",
  function(object){
    standardGeneric("reduction")
  }
)

setGeneric(
  "multiplyOnKoef",
  function(object,koef){
    standardGeneric("multiplyOnKoef")
  }
)

setGeneric(
  "multiplyPolynom",
  function(object){
    standardGeneric("multiplyPolynom")
  }
)

setGeneric(
  "sumEqualsTerms",
  function(object){
    standardGeneric("sumEqualsTerms")
  }
)

setGeneric(
  "getLagranjStringSummand",
  function(object,a){
    standardGeneric("getLagranjStringSummand")
  }
)

setGeneric(
  "lagranj",
  function(object){
    standardGeneric("lagranj")
  }
)

setMethod(  # find quotient for one multiplier from basic polynom
  "symbolQuotient",
  signature = "Interpolation",
  function(object,j,i){
    denominator = 1/(object@x[j] - object@x[i])
    numerator = object@x[i]
    res = c(numerator,denominator)
    return(res)
  }
)


setMethod( #multiply TWO different polynom
  "multiplyPolynom",
  signature = "Interpolation",
  function(object){
    lengthOfL = length(object@l)
    lengthOfTL = length(object@tmpList)
    index=1
    tmpAmount = lengthOfL*lengthOfTL
    newList = list()
    for(i in 1:tmpAmount){
      newList[[i]]=c(0,0)
    }
    for(i in 1:lengthOfL){
      for(j in 1:lengthOfTL){
        newList[[index]][1]=object@l[[i]][1]*object@tmpList[[j]][1]
        newList[[index]][2]=object@l[[i]][2]+object@tmpList[[j]][2]
        index=index+1
      }
    }
    return(newList)
  }
)


setMethod(
  "findBasicPolynom",
  signature = "Interpolation",
  function(object,j){
    if(j!=1){
      x = symbolQuotient(object,j,1)
      q = 1
    } else {
      x = symbolQuotient(object,j,2)
      q = 2
    }
    denominator = x[2]
    lengthOfX = length(object@x)
    for(i in 1:2){
      object@tmpList[[i]]=c(0,0)
    }
    for(i in 1:lengthOfX){
      res=c(0,0)
    }
    object@l[[1]]=c(1,1)
    object@l[[2]]=c(-x[1],0)
    
    for(i in 1:lengthOfX){
      if((i!=j)&&(i!=q)){
        x = symbolQuotient(object,j,i)
        denominator = denominator*x[2]
        object@tmpList[[1]]=c(1,1);
        object@tmpList[[2]]=c(-x[1],0);
        res=multiplyPolynom(object)
        object@l=res
      }
    }
    res[[length(res)+1]]=denominator
    return(res)
  }
)

setMethod( # reduction all zero
  "reduction",
  signature = "Interpolation",
  function(object){
    lengthL= length(object@l)
    i = 1
    while(i<=length(object@l)){
      if(object@l[[i]][1]==0){
        object@l[[i]] = NULL
        lengthL = length(object@l)
      } else  i = i+1
    }
    res = object@l
    return(res)
  }
)

setMethod(# sum equals terms
  "sumEqualsTerms",
  signature = "Interpolation",
  function(object){
    lengthL= length(object@l)
    i = 0
    while(i<lengthL){
      i=i+1
      j=i+1
      while(j<=length(object@l)){
        if(object@l[[i]][2]==object@l[[j]][2]){
          object@l[[i]][1] =  object@l[[i]][1] + object@l[[j]][1]
          object@l[[j]] = NULL
          lengthL = lengthL - 1
        }else j = j + 1
      }
    }
    res = object@l
    return(res)
  }
)

setMethod(
  "findLagranjPolynom",
  signature = "Interpolation",
  function(object){
    res = list()
    for(i in 1:length(object@y)){
      object@l=list()
      tmp = findBasicPolynom(object,i)
      denominator=tmp[[length(tmp)]]
      tmp[[length(tmp)]] = NULL
      object@l = tmp
      object@l = sumEqualsTerms(object)
      object@l = reduction(object)
      denominator = denominator*object@y[i]
      object@l = multiplyOnKoef(object,denominator)
      tmpList = object@l;
      a = 1;
      for(j in (length(res)+1):(length(res)+length(tmpList))){
        res[[j]] = tmpList[[a]]
        a=a+1
      }
    }
    return(res)
  }
)

setMethod(
  "multiplyOnKoef",
  signature = "Interpolation",
  function(object,koef){
    for(i in 1:length(object@l)){
      object@l[[i]][1] = object@l[[i]][1]*koef
    }
    res = object@l
    return(res)
  }
)

setMethod(
  "getLagranjPolynom",
  signature = "Interpolation",
  function(object){
    object@l = findLagranjPolynom(object)
    object@l = sumEqualsTerms(object)
    res = reduction(object)
    return(res)
  }
)

setMethod(
  "getLagranjStringSummand",
  signature = "Interpolation",
  function(object,a){
    res = list()
    for(i in 1:length(a)){
      res[[i]] =  paste("(",paste(a[[i]][1],paste("* x ^",paste(as.numeric(a[[i]][2]),")"))))
    }
    return(res)
  }
)

setMethod(
  "getLagranjStringPolynom",
  signature = "Interpolation",
  function(object,a){
    res = ""
    for(i in 1:length(a)){
      if(i==1) res = paste(res,a[[i]])
      else res = paste(res,paste("+",a[[i]]))
    }
    return(res)
  }
)

setMethod (
  "lagranj",
  signature = "Interpolation",
  function(object){
    listLagranj = getLagranjPolynom(object)
    listLagranj = getLagranjStringSummand(object,listLagranj)
    listLagranj = getLagranjStringPolynom(object,listLagranj)
    return(listLagranj)
  }
)

test = new ("Interpolation",
            y=c(5,0,3,4),
            x=c(-5,-2,1,4)
)


print(lagranj(test))

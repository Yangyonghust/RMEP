#' Multiple imputation methods for multivariate missing data
#'
#' Deal with the missing data of the inputs for the maximum entropy production(MEP) model.
#' The function creates multiple imputations (replacement values) for multivariate missing data.
#' The imputation methods are established based on the mice package:Generates Multivariate Imputations by Chained Equations (MICE).
#' There are mainly six univariate imputation methods available including:Predictive mean matching("pmm"), Weighted predictive mean matching("midastouch"),
#' Classification and regression trees ("cart"), Unconditional mean imputation ("mean"), Imputation of quadratic terms("quadratic") and Linear regression, predicted values ("norm.predict").
#' Users can select method for imputations and the default is Predictive mean matching("pmm").
#'
#' @param data A dataframe of variables but contains missing data for MEP inputs
#' @param method Methods for imputation, including six methods selected
#' @return A dataframe with no missing values
#' @importFrom mice mice complete
#' @examples
#' RMEP_mice(airquality,method= 'pmm')
#' @export
RMEP_mice=function(data,method= 'pmm'){{
  #library(mice)
  completeData=mice(data, m=5, maxit = 50, seed = 500,printFlag = FALSE)  ####  Predictive mean matching
  if(method=="midastouch"){                                   #### Weighted predictive mean matching
    completeData=mice(data, m=5, maxit = 50, method = 'midastouch', seed = 500,printFlag = FALSE)
  }

  else if(method=="cart"){
    completeData=mice(data, m=5, maxit = 50, method = 'cart', seed = 500,printFlag = FALSE)
  }
  else if(method=="mean"){
    completeData=mice(data, m=5, maxit = 50, method = 'mean', seed = 500,printFlag = FALSE)
  }
  else if(method=="quadratic"){
    completeData=mice(data, m=5, maxit = 50, method = 'quadratic', seed = 500,printFlag = FALSE)
  }

  else{
    completeData=mice(data, m=5, maxit = 50, method = 'norm.predict', seed = 500,printFlag = FALSE)
  }
  completeData <- complete(completeData,1)
  return(completeData)
}
}

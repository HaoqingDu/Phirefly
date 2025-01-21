init.prior <- function() {
  
  prior_list <<- list(
    pars = character(),
    dist = character()
  )
}

check.prior <- function() {
  cmd <- paste0(input$pars[i], "<<-", input$dist[i])

  print(cmd)
}

add.prior <- function(parameter, distribution) {
  
  if (! length(parameter) == length(distribution)) {
    stop("Error in prior list: each parameter should have one prior distribution.")
  }
  
  prior_list$pars <<- c(prior_list$pars, parameter)
  prior_list$dist <<- c(prior_list$dist, distribution)
  
  return(
    paste0("Current number of priors: ", length(prior_list$pars))
  )
}


draw.from.prior <- function(input = prior_list) {
  print("Priors:")
  
  if (length) {
    
  }
  
  for (i in 1:length(input$pars)) {
    cmd <- paste0(input$pars[i], "<<-", input$dist[i])
    
    eval(parse(text = cmd))
    
    print(cmd)
  }  
}

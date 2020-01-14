server <- function(input, output, session) {
  
  #####################################################################
  #                   Hist 1
  #####################################################################
  
  output$hist1 <- renderPlot({
    par(mfrow=c(1,2))
    x    <- seq(0,input$N,1)
    y <- dbinom(x, input$N, input$nullvalue)
    
    #decisionthreshold pvalue
    ThresholdRight = (qbinom(input$nullvalue* input$critp,input$N,input$nullvalue, lower.tail = FALSE))+1 # inverse cumulative probability for p=0.05, what should be the minimal number of heads to conclude significance?
    ThresholdLeft = (qbinom(input$nullvalue* input$critp,input$N,input$nullvalue, lower.tail = TRUE))-1 # inverse cumulative probability for p=0.05, what should be the minimal number of heads to conclude significance?
    RejL = 1+ ThresholdLeft #how many y's are in the left rejection area
    RejR = input$N - ThresholdRight #how many y's are in the right rejection
    
    ### Make the prop false alarm red & power light darkgrey in plot
    RejRVec = seq(ThresholdRight,input$N,by=1) # a vector containing all the number of 'heads' that would be significant
    RejLVec = seq(0,ThresholdLeft,by=1) # a vector containing all the number of 'heads' that would be significant
    
    barplot(y, col = c((rep((rgb(1, 0, 0,0.5)), RejL)), 
                       (rep("darkgray", (input$N-RejL-RejR))), 
                       (rep((rgb(1, 0, 0,0.5)), RejR))), 
            border = 'grey',
            main = "What outcome would be significant?",
            ylab = "Probability Density", 
            xlab = "Outcome - number of heads",
            names.arg = x)
    
    ynull <- dbinom(x, input$N, input$nullvalue) #nullhypothesis
    yalt <- (dbeta((x/input$N),input$aprior,input$bprior))/(input$N+1) #althypothesis
    
    if (input$ratio == "BF10") {
      vecratioaltnull=yalt/ynull #calculate the BF for every prop
      ratioaltnullpvalueR=yalt[ThresholdRight+1]/ynull[ThresholdRight+1] #what is the critical BF that corresponds to same cut off point as the p-value?
      ratioaltnullpvalueL=yalt[ThresholdLeft+1]/ynull[ThresholdLeft+1] #what is the critical BF that corresponds to same cut off point as the p-value?
    }
    
    if (input$ratio == "BF01") {
      vecratioaltnull=ynull/yalt #calculate the BF for every prop
      ratioaltnullpvalueR=ynull[ThresholdRight+1]/yalt[ThresholdRight+1] #what is the critical BF that corresponds to same cut off point as the p-value?
      ratioaltnullpvalueL=ynull[ThresholdLeft+1]/yalt[ThresholdLeft+1]#what is the critical BF that corresponds to same cut off point as the p-value?
    }
    
    barplot(vecratioaltnull, ylim=c(0,input$ymax),xpd=F, col = c((rep((rgb(1, 0, 0,0.5)), RejL)), 
                                                                 (rep("darkgray", (input$N-(RejL+RejR)))), 
                                                                 (rep((rgb(1, 0, 0,0.5)), RejR))), 
            border = 'grey',
            main = "When would we conclude there's an effect?",
            ylab = "BayesFactor",
            xlab = "Outcome - number of heads",
            names.arg = x)
  })
  
  #######################################################################################
  #Give corresponding BF
  #######################################################################################
  
  output$summary <- renderText({
    
    # if 1
    
    if (input$aprior == input$bprior) {
      
      x=qbinom(0.5*input$critp, size=input$N, prob = input$nullvalue) #from alfa to minimal observed number of heads to be significant
      if (2*(pbinom(x, input$N, input$nullvalue)) > input$critp) {x= qbinom(0.5*input$critp, size=input$N, prob = input$nullvalue)-1}
      
      if (input$ratio=="BF01") {bfcrit= 1/(( exp( lbeta(input$aprior+x,input$bprior+input$N-x) - lbeta(input$aprior,input$bprior) ) / ( input$nullvalue^x * (1.0-input$nullvalue)^(input$N-x) ) ))} 
      if (input$ratio=="BF10") {bfcrit= ( exp( lbeta(input$aprior+x,input$bprior+input$N-x) - lbeta(input$aprior,input$bprior) ) / ( input$nullvalue^x * (1.0-input$nullvalue)^(input$N-x) ) )   }
      
      paste("The corresponding Bayes Factor is", round(bfcrit,2))}
    
    # if 2.1
    
    else { 
      xl=qbinom(0.5*input$critp, size=input$N, prob = input$nullvalue, lower.tail=TRUE) #from alfa to maximal observed number of heads to be significant (left tail)
      if (2*(pbinom(xl, input$N, input$nullvalue)) > input$critp) {xl= qbinom(0.5*input$critp, size=input$N, prob = input$nullvalue, lower.tail=TRUE)-1}
      
      xu=qbinom(0.5*input$critp, size=input$N, prob = input$nullvalue, lower.tail=FALSE) #from alfa to minimal observed number of heads to be significant (right tail)
      if (2*(pbinom(xu, input$N, input$nullvalue)) > input$critp) {xu= qbinom(0.5*input$critp, size=input$N, prob = input$nullvalue, lower.tail=FALSE)+1}
      
      if (input$ratio=="BF10") {
        A = ( exp( lbeta(input$aprior+xl,input$bprior+input$N-xl) - lbeta(input$aprior,input$bprior) ) / ( input$nullvalue^xl * (1.0-input$nullvalue)^(input$N-xl) ) ) 
        lA = ( exp( lbeta(input$aprior+(xl-1),input$bprior+input$N-(xl-1)) - lbeta(input$aprior,input$bprior) ) / ( input$nullvalue^(xl-1) * (1.0-input$nullvalue)^(input$N-(xl-1)) ) ) 
        
        B = ( exp( lbeta(input$aprior+xu,input$bprior+input$N-xu) - lbeta(input$aprior,input$bprior) ) / ( input$nullvalue^xu * (1.0-input$nullvalue)^(input$N-xu) ) ) 
        uB = ( exp( lbeta(input$aprior+(xu+1),input$bprior+input$N-(xu+1)) - lbeta(input$aprior,input$bprior) ) / ( input$nullvalue^(xu+1) * (1.0-input$nullvalue)^(input$N-(xu+1)) ) ) 
        
        
        if (A>lA) {signa = "<"}
        if (A<lA) {signa = ">"}
        
        if (B<uB) {signb = ">"}
        if (B>uB) {signb = "<"}
        
        bfcrit <- (paste("The Bayes Factor graph is not symmetrical, to obtain the same rejection regions we need two different Bayes Factors: BFleft", signa, round(A,2), "and BFright", signb, round(B,2))) }
      
      # if 2.2
      
      if (input$ratio=="BF01") {
        A = (1/( exp( lbeta(input$aprior+xl,input$bprior+input$N-xl) - lbeta(input$aprior,input$bprior) ) / ( input$nullvalue^xl * (1.0-input$nullvalue)^(input$N-xl) ) ) )
        lA = (1/( exp( lbeta(input$aprior+(xl-1),input$bprior+input$N-(xl-1)) - lbeta(input$aprior,input$bprior) ) / ( input$nullvalue^(xl-1) * (1.0-input$nullvalue)^(input$N-(xl-1)) ) ) )
        
        B = (1/( exp( lbeta(input$aprior+xu,input$bprior+input$N-xu) - lbeta(input$aprior,input$bprior) ) / ( input$nullvalue^xu * (1.0-input$nullvalue)^(input$N-xu) ) ) )
        uB = (1/( exp( lbeta(input$aprior+(xu+1),input$bprior+input$N-(xu+1)) - lbeta(input$aprior,input$bprior) ) / ( input$nullvalue^(xu+1) * (1.0-input$nullvalue)^(input$N-(xu+1)) ) ) )
        
        
        if (A>lA) {signa = "<"}
        if (A<lA) {signa = ">"}
        
        if (B<uB) {signb = ">"}
        if (B>uB) {signb = "<"}
        
        bfcrit <- (paste("The Bayes Factor graph is not symmetrical, to obtain the same rejection regions we need two different Bayes Factors: BFleft", signa, round(A, 2), "and BFright", signb, round(B,2)))  #### output optie 2
      }
      bfcrit} 
    
  })
}
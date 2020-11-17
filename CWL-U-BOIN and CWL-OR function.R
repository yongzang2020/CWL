
library(nleqslv)
library(BOIN)
library(rjags)

########### The following part is for CWL U-BOIN function #############

######## We use following function to get the suitable parameter of weibull distribution for U-BOIN simualtion
#' @param U_T is the assessment period for DLT
#' @param v_T is the interim period for DLT
#' @param U_E is the assessment period for efficacy
#' @param v_E is the interim period for efficacy
#' @param DLT is DLT rate curve for each dose level within assessment period
#' @param eff is efficacy rate curve for each dose level within assessment period
#' @param event_rate is the efficacy and DLT event rate happens within V_T & V_E
#' @return the suitable parameters of weibull distribution for DLT and efficacy 

##############
para_finding_BOIN<-function(U_T=3,V_T=2,U_E=3,V_E=2, event_rate=0.6, DLT=DLT,eff=efficacy){
#######
para_U<- matrix(0,nrow=5,ncol=6)
colnames(para_U)<-c("p_k","p_lambda","q0_k","q0_lambda","q1_k","q1_lambda")

for (i in 1:length(DLT)) {
  
  ###### DLT rate
  p_U<-DLT[i]
  p_Z<-DLT[i]*event_rate
  
  ###### Efficacy
  q0_U<-eff[i]
  #######
  q0_Z<-q0_U*event_rate
  
  #### pars[1]=shape=k pars[2]=lambda=scale
  f_1<-function(pars){
    
    y<- numeric(2)
    y[1]<-pweibull(U_T, pars[1],pars[2], lower.tail = TRUE, log.p = FALSE) - p_U
    y[2]<-pweibull(V_T, pars[1],pars[2], lower.tail = TRUE, log.p = FALSE) - p_Z
    
    return (y)
  }
  
  
  ### xstart is start point , nleqslv is the function to find parameter
  xstart <- c(2,3)
  
  outcome1<-nleqslv(xstart, fn=f_1, method="Newton",global="qline",control = list(allowSingular=TRUE))
  
  
  ### get the result for toxiciy's parameters k and lambda
  k<-outcome1$x[1]
  lambda<-outcome1$x[2]
  
  
  #### pars[1]=shape=k_q0 pars[2]=scale=lambda_q0
  f_2<-function(pars){
    
    y<- numeric(2)
    y[1]<-pweibull(U_E, pars[1],pars[2], lower.tail = TRUE, log.p = FALSE) - q0_U
    y[2]<-pweibull(V_E, pars[1],pars[2], lower.tail = TRUE, log.p = FALSE) - q0_Z
    
    return (y)
  }
  
  
  xstart<-c(2,3)
  
  outcome2<-nleqslv(xstart, fn=f_2, method="Newton",global="qline",control = list(allowSingular=TRUE))
 
  ### get the result for conditional efficacy's parameters k and lambda
  k_q0<-outcome2$x[1]
  lambda_q0<-outcome2$x[2]
  
  k_q1<-k_q0
  lambda_q1<-lambda_q0
  
  ###### Save the para
  para_U[i,]<-c(k,lambda,k_q0,lambda_q0,k_q1,lambda_q1)
  
}

return(para_U)
}

######### We use folloing function to do the simulation for CWL U-BOIN method
#' @param para is the initial parameter result geeting from above function 
#' @param U_T is the assessment period for DLT
#' @param v_T is the interim period for DLT
#' @param U_E is the assessment period for efficacy
#' @param v_E is the interim period for efficacy
#' @param pi_T is pre-setting DLT rate upper boundary for the medicine
#' @param pi_E is pre-setting efficacy rate lower boundary for the medicine
#' @param p_cut is the posterior distribution cut-off for DLT rate
#' @param q_cut is the posterior distribution cut-off for efficacy rate
#' @param no_cohort is the totol cohorts involved in one clinicl trial
#' @param cohort_size is the cohort size for one cohort
#' @param score is the mark score (omega) for four different status of patients
#' @param M is the total simualtion times.

CWL_U_BOIN_OC<-function(para_U=para,phi_0=0.3,U_T=3,U_E=3,V_T=2,V_E=2,p_cut=0.95,q_cut=0.9,
                        pi_T=0.3,pi_E=0.2, no_cohort=20, cohort_size=3,score=c(0,30,50,100), M=1000) {

##### Set raw parameters

phi<-(phi_0-0.05)

phi_1<-0.6*phi
phi_2<-1.4*phi

##### Calculate the boundary according to phis
lambda_1<-log((1-phi_1)/(1-phi))/log(phi*(1-phi_1)/phi_1/(1-phi))

lambda_2<-log((1-phi)/(1-phi_2))/log(phi_2*(1-phi)/phi/(1-phi_2))

#### Prior is Beta(0.5,0.5)

######## Stage II function

#### Write the model with 5 dose levels
model_U<- "model{

#   # Likelihood
y1~dmulti(R1,n1)
y2~dmulti(R2,n2)
y3~dmulti(R3,n3)
y4~dmulti(R4,n4)
y5~dmulti(R5,n5)
y6~dmulti(R6,n6)
y7~dmulti(R7,n7)
y8~dmulti(R8,n8)
y9~dmulti(R9,n9)
y10~dmulti(R10,n10)



#   # Prior
#   # 1 is toxicity and efficacy; 2 is toxicity and non-efficacy
#   # 3 is non-toxiciy and efficacy;4 is non-toxicity and non-efficacy;
##### Fully observed
R1[1]<- p_1*q_11
R1[2]<- p_1*(1-q_11)
R1[3]<- (1-p_1)*q_01
R1[4]<- (1-p_1)*(1-q_01)


R2[1]<- p_2*q_12
R2[2]<- p_2*(1-q_12)
R2[3]<- (1-p_2)*q_02
R2[4]<- (1-p_2)*(1-q_02)


R3[1]<- p_3*q_13
R3[2]<- p_3*(1-q_13)
R3[3]<- (1-p_3)*q_03
R3[4]<- (1-p_3)*(1-q_03)


R4[1]<- p_4*q_14
R4[2]<- p_4*(1-q_14)
R4[3]<- (1-p_4)*q_04
R4[4]<- (1-p_4)*(1-q_04)


R5[1]<- p_5*q_15
R5[2]<- p_5*(1-q_15)
R5[3]<- (1-p_5)*q_05
R5[4]<- (1-p_5)*(1-q_05)


##### Partly missing----2 tau
R6[1]<- (V_E/U_E)*q_11*(V_T/U_T)*p_1
R6[2]<- (V_T/U_T)*p_1*(1-(V_E/U_E)*q_11)
R6[3]<- (1-(V_T/U_T)*p_1)*(V_E/U_E)*(q_11*((1-V_T/U_T)*p_1/((1-V_T/U_T)*p_1+1-p_1))+q_01*((1-p_1)/((1-V_T/U_T)*p_1+1-p_1)))
R6[4]<- (1-(V_T/U_T)*p_1)*(1-(V_E/U_E)*(q_11*((1-V_T/U_T)*p_1/((1-V_T/U_T)*p_1+1-p_1))+q_01*((1-p_1)/((1-V_T/U_T)*p_1+1-p_1))))


R7[1]<- (V_E/U_E)*q_12*(V_T/U_T)*p_2
R7[2]<- (V_T/U_T)*p_2*(1-(V_E/U_E)*q_12)
R7[3]<- (1-(V_T/U_T)*p_2)*(V_E/U_E)*(q_12*((1-V_T/U_T)*p_2/((1-V_T/U_T)*p_2+1-p_2))+q_02*((1-p_2)/((1-V_T/U_T)*p_2+1-p_2)))
R7[4]<- (1-(V_T/U_T)*p_2)*(1-(V_E/U_E)*(q_12*((1-V_T/U_T)*p_2/((1-V_T/U_T)*p_2+1-p_2))+q_02*((1-p_2)/((1-V_T/U_T)*p_2+1-p_2))))


R8[1]<- (V_E/U_E)*q_13*(V_T/U_T)*p_3
R8[2]<- (V_T/U_T)*p_3*(1-(V_E/U_E)*q_13)
R8[3]<- (1-(V_T/U_T)*p_3)*(V_E/U_E)*(q_13*((1-V_T/U_T)*p_3/((1-V_T/U_T)*p_3+1-p_3))+q_03*((1-p_3)/((1-V_T/U_T)*p_3+1-p_3)))
R8[4]<- (1-(V_T/U_T)*p_3)*(1-(V_E/U_E)*(q_13*((1-V_T/U_T)*p_3/((1-V_T/U_T)*p_3+1-p_3))+q_03*((1-p_3)/((1-V_T/U_T)*p_3+1-p_3))))


R9[1]<- (V_E/U_E)*q_14*(V_T/U_T)*p_4
R9[2]<- (V_T/U_T)*p_4*(1-(V_E/U_E)*q_14)
R9[3]<- (1-(V_T/U_T)*p_4)*(V_E/U_E)*(q_14*((1-V_T/U_T)*p_4/((1-V_T/U_T)*p_4+1-p_4))+q_04*((1-p_4)/((1-V_T/U_T)*p_4+1-p_4)))
R9[4]<- (1-(V_T/U_T)*p_4)*(1-(V_E/U_E)*(q_14*((1-V_T/U_T)*p_4/((1-V_T/U_T)*p_4+1-p_4))+q_04*((1-p_4)/((1-V_T/U_T)*p_4+1-p_4))))


R10[1]<- (V_E/U_E)*q_15*(V_T/U_T)*p_5
R10[2]<- (V_T/U_T)*p_5*(1-(V_E/U_E)*q_15)
R10[3]<- (1-(V_T/U_T)*p_5)*(V_E/U_E)*(q_15*((1-V_T/U_T)*p_5/((1-V_T/U_T)*p_5+1-p_5))+q_05*((1-p_5)/((1-V_T/U_T)*p_5+1-p_5)))
R10[4]<- (1-(V_T/U_T)*p_5)*(1-(V_E/U_E)*(q_15*((1-V_T/U_T)*p_5/((1-V_T/U_T)*p_5+1-p_5))+q_05*((1-p_5)/((1-V_T/U_T)*p_5+1-p_5))))

######## Toxicity and efficacy planned follow-up time and fixed follow-up time

U_E<-3
U_T<-3

####### actual follow-up time
V_E<-2
V_T<-2


########
### Dose 1
p_1~dbeta(1,1)
q_01~dbeta(1,1)
q_11~dbeta(1,1)


### Dose 2
p_2~dbeta(1,1)
q_02~dbeta(1,1)
q_12~dbeta(1,1)


### Dose 3
p_3~dbeta(1,1)
q_03~dbeta(1,1)
q_13~dbeta(1,1)

### Dose 4
p_4~dbeta(1,1)
q_04~dbeta(1,1)
q_14~dbeta(1,1)

### Dose 5
p_5~dbeta(1,1)
q_05~dbeta(1,1)
q_15~dbeta(1,1)

}"



##### Original safety rule parameters which should be adjusted 
##### Only used in stage II

options(digits=5)
####Big simulation

#### Create the Big summary table
Big_summary<-matrix(0,M,12)

dose_history_sum<-matrix(0,M,no_cohort+1)


for (k in 1:M) {
  
  ### create the accrue data frame, accrue_s_a, accrue_delta_a
  accrue_D<-matrix(8,1,7)
  accrue_D<-accrue_D[-1,]
  
  #### Rest the parameter
  dose<-1
  
  dose_history<-c(1,rep(0,no_cohort))
  
  stage<-rep(0,no_cohort)
  
  summary_dose<-rep(0,5)
  
  time<-0
  
  early_stopping<- 0
  
  weight<-V_T/U_T
  
  admissible<-rep(1,5)
  
  for (l in 1:no_cohort) {
    
    #### Tempory data
    temp_D<-matrix(0,cohort_size,7)
    
    temp_D[,1]<-rep(dose,cohort_size)
    
    for (i in 1:cohort_size){
      
      temp_D[i,2]<-rweibull(1,shape = para_U[dose,1], scale = para_U[dose,2])
      temp_D[i,3]<-(temp_D[i,2]<=V_T)
      temp_D[i,4]<-(temp_D[i,2]<=U_T)
      
      #### Based on toxicity to simulate efficacy  
      temp_D[i,5]<-ifelse((temp_D[i,2]>=U_T),rweibull(1,shape = para_U[dose,3], scale = para_U[dose,4]),
                          rweibull(1,shape = para_U[dose,5], scale = para_U[dose,6]))
      
      temp_D[i,6]<-(temp_D[i,5]<=V_E)
      temp_D[i,7]<-(temp_D[i,5]<=U_E)
    }
    
    
    ####accruetive data
    accrue_D<-rbind(accrue_D,temp_D)
    accrue_fake<- accrue_D[1:(nrow(accrue_D)-cohort_size),]
    
    ###### Judgement function for identifying stage I and stage II #####
    
    for (c in 1:5) {    
      
      summary_dose[c]<-sum(accrue_D[,1]==c)
      
    }
    ##### if any(summary_dose>=12) then in stage II else, keep stage I
    if (max(summary_dose)<12) 
    { 
      ##### Track stage
      stage[l]<-1
      
      #### Calculate toxicity reponse and total patients number
      {
        if (l==1)
        {n_toxi<-c(sum(temp_D[,3]),0,0,0,0)
        
        n_m_toxi<-c(sum(temp_D[,3]==0)*weight,0,0,0,0)
        }
        else if ((l>1)*(l<no_cohort))
        {n_toxi<- c(sum((accrue_fake[,1]==1)*(accrue_fake[,4]))+sum((temp_D[,3]==1)*(temp_D[,1]==1)),
                    sum((accrue_fake[,1]==2)*(accrue_fake[,4]))+sum((temp_D[,3]==1)*(temp_D[,1]==2)),
                    sum((accrue_fake[,1]==3)*(accrue_fake[,4]))+sum((temp_D[,3]==1)*(temp_D[,1]==3)),
                    sum((accrue_fake[,1]==4)*(accrue_fake[,4]))+sum((temp_D[,3]==1)*(temp_D[,1]==4)),
                    sum((accrue_fake[,1]==5)*(accrue_fake[,4]))+sum((temp_D[,3]==1)*(temp_D[,1]==5)))
        
        n_m_toxi<-c(sum((accrue_fake[,1]==1)*(accrue_fake[,4]==0))+sum((temp_D[,3]==0)*(temp_D[,1]==1))*weight,
                    sum((accrue_fake[,1]==2)*(accrue_fake[,4]==0))+sum((temp_D[,3]==0)*(temp_D[,1]==2))*weight,
                    sum((accrue_fake[,1]==3)*(accrue_fake[,4]==0))+sum((temp_D[,3]==0)*(temp_D[,1]==3))*weight,
                    sum((accrue_fake[,1]==4)*(accrue_fake[,4]==0))+sum((temp_D[,3]==0)*(temp_D[,1]==4))*weight,
                    sum((accrue_fake[,1]==5)*(accrue_fake[,4]==0))+sum((temp_D[,3]==0)*(temp_D[,1]==5))*weight)
        }
        else if (l==no_cohort)
        {n_toxi<- c(sum((accrue_D[,1]==1)*(accrue_D[,4])),
                    sum((accrue_D[,1]==2)*(accrue_D[,4])),
                    sum((accrue_D[,1]==3)*(accrue_D[,4])),
                    sum((accrue_D[,1]==4)*(accrue_D[,4])),
                    sum((accrue_D[,1]==5)*(accrue_D[,4])))
        n_m_toxi<-c(sum((accrue_D[,1]==1)*(accrue_D[,4]==0)),
                    sum((accrue_D[,1]==2)*(accrue_D[,4]==0)),
                    sum((accrue_D[,1]==3)*(accrue_D[,4]==0)),
                    sum((accrue_D[,1]==4)*(accrue_D[,4]==0)),
                    sum((accrue_D[,1]==5)*(accrue_D[,4]==0)))
        }
        
        n_dose <- n_toxi+ n_m_toxi
        p_estimat<-n_toxi/n_dose
      }
      #### Get the a accumulative time
      { if (l<no_cohort)
        time<-time+V_T
        else if ((l=no_cohort)*(sum(temp_D[,4])==cohort_size))
          time<-time+max(temp_D[,2])
        else if ((l=no_cohort)*(sum(temp_D[,4])<cohort_size))
          time<-time+U_T
      }
      
      #### Setting Early stopping Rule
      {if ((l<no_cohort)*(dose==1)*(pbeta((pi_T-0.05),0.5+n_toxi[dose],0.5+(n_dose[dose]-n_toxi[dose]),lower.tail = FALSE)>p_cut))
        
      {early_stopping <- 1
      
      dose<-0
      
      break}
        #### Get current admissible set & break condition
        else if ((l==no_cohort)*(dose==1)*(pbeta((pi_T-0.05),0.5+n_toxi[dose],0.5+(n_dose[dose]-n_toxi[dose]),lower.tail = FALSE)>p_cut))
          
        {early_stopping <-0
        
        dose<-0
        
        break}
        
        else if (pbeta((pi_T-0.05),0.5+n_toxi[dose],0.5+(n_dose[dose]-n_toxi[dose]),lower.tail = FALSE)>p_cut)
          
          admissible<-c(rep(1,(dose-1)),rep(0,(6-dose)))
      }
      
      #### Get maximum dose in admissible set
      max_admissible<-max(admissible*c(1,2,3,4,5))
      
      #### Next dose selection
      
      {if ((p_estimat[dose]>=lambda_2)*(dose==1))
        dose<-1
        else if (p_estimat[dose]>=lambda_2)
          dose<-dose-1
        else if ((p_estimat[dose]<=lambda_1)*(dose==max_admissible))
          dose<-dose
        else if (p_estimat[dose]<=lambda_1)
          dose<-min(dose+1,max_admissible)
        else dose<-dose}
    }
    
    ############ Another judge statement for stage II ###########    
    else if (any(summary_dose>=12))
      
    {##### Track stage
      stage[l]<-2
      y<-matrix(0,10,4)
      
      if (l<no_cohort)
      { for (j in 1:5){
        y[j,1]<-sum(accrue_D[1:(3*l-3),4]*accrue_D[1:(3*l-3),7]*(accrue_D[1:(3*l-3),1]==j))
        y[j,2]<-sum(accrue_D[1:(3*l-3),4]*(1-accrue_D[1:(3*l-3),7])*(accrue_D[1:(3*l-3),1]==j))
        y[j,3]<-sum((1-accrue_D[1:(3*l-3),4])*accrue_D[1:(3*l-3),7]*(accrue_D[1:(3*l-3),1]==j))
        y[j,4]<-sum((1-accrue_D[1:(3*l-3),4])*(1-accrue_D[1:(3*l-3),7])*(accrue_D[1:(3*l-3),1]==j))
      }
        for (j in 6:10){
          y[j,1]<-sum(accrue_D[(3*l-2):(3*l),3]*accrue_D[(3*l-2):(3*l),6]*(accrue_D[(3*l-2):(3*l),1]==(j-5)))
          y[j,2]<-sum(accrue_D[(3*l-2):(3*l),3]*(1-accrue_D[(3*l-2):(3*l),6])*(accrue_D[(3*l-2):(3*l),1]==(j-5)))
          y[j,3]<-sum((1-accrue_D[(3*l-2):(3*l),3])*accrue_D[(3*l-2):(3*l),6]*(accrue_D[(3*l-2):(3*l),1]==(j-5)))
          y[j,4]<-sum((1-accrue_D[(3*l-2):(3*l),3])*(1-accrue_D[(3*l-2):(3*l),6])*(accrue_D[(3*l-2):(3*l),1]==(j-5)))
        }
      }
      ######
      else if (l==no_cohort)
      {for (j in 1:5){
        y[j,1]<-sum(accrue_D[,4]*accrue_D[,7]*(accrue_D[,1]==j))
        y[j,2]<-sum(accrue_D[,4]*(1-accrue_D[,7])*(accrue_D[,1]==j))
        y[j,3]<-sum((1-accrue_D[,4])*accrue_D[,7]*(accrue_D[,1]==j))
        y[j,4]<-sum((1-accrue_D[,4])*(1-accrue_D[,7])*(accrue_D[,1]==j))
      }
      }
      
      y1<-y[1,]
      y2<-y[2,]
      y3<-y[3,]
      y4<-y[4,]
      y5<-y[5,]
      y6<-y[6,]
      y7<-y[7,]
      y8<-y[8,]
      y9<-y[9,]
      y10<-y[10,]
      #####
      n1<-sum(y1)
      n2<-sum(y2)
      n3<-sum(y3)
      n4<-sum(y4)
      n5<-sum(y5)
      #####
      n6<-sum(y6)
      n7<-sum(y7)
      n8<-sum(y8)
      n9<-sum(y9)
      n10<-sum(y10)
      
      ######
      model <- jags.model(textConnection(model_U), data = list(y1=y1,y2=y2,y3=y3,y4=y4,y5=y5,y6=y6,y7=y7,y8=y8,y9=y9,y10=y10, n1=n1,n2=n2,n3=n3,n4=n4,n5=n5,n6=n6,n7=n7,n8=n8,n9=n9,n10=n10),quiet=T)
      
      update(model, 1000, progress.bar="none")
      
      S<-5000
      
      samp <- jags.samples(model, n.iter=S, variable.names=c("p_1","p_2","p_3","p_4","p_5","q_01","q_02","q_03","q_04","q_05","q_11","q_12","q_13","q_14","q_15"), 
                           progress.bar="none")
      
      p_1<- samp$p_1[1:S]
      p_2<- samp$p_2[1:S]
      p_3<- samp$p_3[1:S]
      p_4<- samp$p_4[1:S]
      p_5<- samp$p_5[1:S]
      
      q_01<- samp$q_01[1:S]
      q_02<- samp$q_02[1:S]
      q_03<- samp$q_03[1:S]
      q_04<- samp$q_04[1:S]
      q_05<- samp$q_05[1:S]
      
      q_11<- samp$q_11[1:S]
      q_12<- samp$q_12[1:S]
      q_13<- samp$q_13[1:S]
      q_14<- samp$q_14[1:S]
      q_15<- samp$q_15[1:S]
      
      
      ##### Combine the original result
      sim_result<-cbind(p_1,q_01,q_11,p_1*(1-q_11),(1-p_1)*(1-q_01),p_1*q_11,(1-p_1)*q_01,
                        p_2,q_02,q_12,p_2*(1-q_12),(1-p_2)*(1-q_02),p_2*q_12,(1-p_2)*q_02,
                        p_3,q_03,q_13,p_3*(1-q_13),(1-p_3)*(1-q_03),p_3*q_13,(1-p_3)*q_03,
                        p_4,q_04,q_14,p_4*(1-q_14),(1-p_4)*(1-q_04),p_4*q_14,(1-p_4)*q_04,
                        p_5,q_05,q_15,p_5*(1-q_15),(1-p_5)*(1-q_05),p_5*q_15,(1-p_5)*q_05)
      
      
      sum_sim_result<-matrix(0,5,4)
      
      for (j in 1:5) {
        
        for (t in 1:4) {
          
          sum_sim_result[j,t]<-mean(sim_result[,(j-1)*7+(t+3)])
          
        }  
      }
      
      ##### update admissible set
      
      { ##### Toxicity & efficacy when l==1
        
        if (l==1)
        {n_toxi<-c(sum(temp_D[,3]),0,0,0,0)
        
        n_m_toxi<-c(sum(temp_D[,3]==0)*weight,0,0,0,0)
        
        n_effi<-c(sum(temp_D[,6]),0,0,0,0)
        
        n_m_effi<-c(sum(temp_D[,6]==0)*weight,0,0,0,0)
        }
        
        else if ((l>1)*(l<no_cohort))
          ##### Toxicity 
        {n_toxi<- c(sum((accrue_fake[,1]==1)*(accrue_fake[,4]))+sum((temp_D[,3]==1)*(temp_D[,1]==1)),
                    sum((accrue_fake[,1]==2)*(accrue_fake[,4]))+sum((temp_D[,3]==1)*(temp_D[,1]==2)),
                    sum((accrue_fake[,1]==3)*(accrue_fake[,4]))+sum((temp_D[,3]==1)*(temp_D[,1]==3)),
                    sum((accrue_fake[,1]==4)*(accrue_fake[,4]))+sum((temp_D[,3]==1)*(temp_D[,1]==4)),
                    sum((accrue_fake[,1]==5)*(accrue_fake[,4]))+sum((temp_D[,3]==1)*(temp_D[,1]==5)))
        
        n_m_toxi<-c(sum((accrue_fake[,1]==1)*(accrue_fake[,4]==0))+sum((temp_D[,3]==0)*(temp_D[,1]==1))*weight,
                    sum((accrue_fake[,1]==2)*(accrue_fake[,4]==0))+sum((temp_D[,3]==0)*(temp_D[,1]==2))*weight,
                    sum((accrue_fake[,1]==3)*(accrue_fake[,4]==0))+sum((temp_D[,3]==0)*(temp_D[,1]==3))*weight,
                    sum((accrue_fake[,1]==4)*(accrue_fake[,4]==0))+sum((temp_D[,3]==0)*(temp_D[,1]==4))*weight,
                    sum((accrue_fake[,1]==5)*(accrue_fake[,4]==0))+sum((temp_D[,3]==0)*(temp_D[,1]==5))*weight)
        
        ##### Efficacy
        n_effi<- c(sum((accrue_fake[,1]==1)*(accrue_fake[,7]))+sum((temp_D[,6]==1)*(temp_D[,1]==1)),
                   sum((accrue_fake[,1]==2)*(accrue_fake[,7]))+sum((temp_D[,6]==1)*(temp_D[,1]==2)),
                   sum((accrue_fake[,1]==3)*(accrue_fake[,7]))+sum((temp_D[,6]==1)*(temp_D[,1]==3)),
                   sum((accrue_fake[,1]==4)*(accrue_fake[,7]))+sum((temp_D[,6]==1)*(temp_D[,1]==4)),
                   sum((accrue_fake[,1]==5)*(accrue_fake[,7]))+sum((temp_D[,6]==1)*(temp_D[,1]==5)))
        
        n_m_effi<-c(sum((accrue_fake[,1]==1)*(accrue_fake[,7]==0))+sum((temp_D[,6]==0)*(temp_D[,1]==1))*weight,
                    sum((accrue_fake[,1]==2)*(accrue_fake[,7]==0))+sum((temp_D[,6]==0)*(temp_D[,1]==2))*weight,
                    sum((accrue_fake[,1]==3)*(accrue_fake[,7]==0))+sum((temp_D[,6]==0)*(temp_D[,1]==3))*weight,
                    sum((accrue_fake[,1]==4)*(accrue_fake[,7]==0))+sum((temp_D[,6]==0)*(temp_D[,1]==4))*weight,
                    sum((accrue_fake[,1]==5)*(accrue_fake[,7]==0))+sum((temp_D[,6]==0)*(temp_D[,1]==5))*weight)
        }
        
        else if (l==no_cohort)
        {n_toxi<- c(sum((accrue_D[,1]==1)*(accrue_D[,4])),
                    sum((accrue_D[,1]==2)*(accrue_D[,4])),
                    sum((accrue_D[,1]==3)*(accrue_D[,4])),
                    sum((accrue_D[,1]==4)*(accrue_D[,4])),
                    sum((accrue_D[,1]==5)*(accrue_D[,4])))
        
        n_m_toxi<-c(sum((accrue_D[,1]==1)*(accrue_D[,4]==0)),
                    sum((accrue_D[,1]==2)*(accrue_D[,4]==0)),
                    sum((accrue_D[,1]==3)*(accrue_D[,4]==0)),
                    sum((accrue_D[,1]==4)*(accrue_D[,4]==0)),
                    sum((accrue_D[,1]==5)*(accrue_D[,4]==0)))
        
        n_effi<- c(sum((accrue_D[,1]==1)*(accrue_D[,7])),
                   sum((accrue_D[,1]==2)*(accrue_D[,7])),
                   sum((accrue_D[,1]==3)*(accrue_D[,7])),
                   sum((accrue_D[,1]==4)*(accrue_D[,7])),
                   sum((accrue_D[,1]==5)*(accrue_D[,7])))
        
        n_m_effi<-c(sum((accrue_D[,1]==1)*(accrue_D[,7]==0)),
                    sum((accrue_D[,1]==2)*(accrue_D[,7]==0)),
                    sum((accrue_D[,1]==3)*(accrue_D[,7]==0)),
                    sum((accrue_D[,1]==4)*(accrue_D[,7]==0)),
                    sum((accrue_D[,1]==5)*(accrue_D[,7]==0)))
        }
        ##### n_dose is the NO of patient of each dose levels
        n_dose <- n_toxi+ n_m_toxi
        
        p_estimat<-n_toxi/n_dose
      }
      
      #### Setting admissible set
      {if (pbeta((pi_T-0.05),0.5+n_toxi[dose],0.5+n_m_toxi[dose],lower.tail = FALSE)>p_cut)
        admissible[dose]<-0
        else if (pbeta(pi_E,0.5+n_effi[dose],0.5+n_m_effi[dose],lower.tail = TRUE)>q_cut)
          admissible[dose]<-0
      }
      
      ###### Calculate the utility
      
      utility<-rep(0,5)
      
      for (i in 1:5) {
        
        utility[i]<- score[1]*sum_sim_result[i,1]+score[2]*sum_sim_result[i,2]+score[3]*sum_sim_result[i,3]+score[4]*sum_sim_result[i,4]
        
      }
      
      ######
      ##### This code without stopping rule and time account
      
      max_dose<-max(dose_history)
      
      max_admissible<-max(admissible*c(1,2,3,4,5))
      
      if (sum(admissible)==0)
      { dose<-0
      
      if(l==no_cohort)
      {early_stopping <- 0
      time<- time+min(max(temp_D[,2],temp_D[,5]),max(U_T,U_E))
      }
      else if (l<no_cohort) {
        early_stopping <- 1
        time<- time+min(max(temp_D[,2],temp_D[,5]),max(V_T,V_E))
      }
      break}
      
      ###### Must change: Force escalate dose to dose+1 phi_2
      else if ((max_dose<5)*(p_estimat[max_dose]<phi_2))
        
        dose<-min(max_dose+1,5)
      
      else if (sum(admissible)>0)
      { if (which.max(admissible*utility)>(max_dose+1))
        dose<-min(max_dose+1,5)
      else if (which.max(admissible*utility)<=(max_dose+1))
        dose<-which.max(admissible*utility)
      }
      
      ######################
      ##### Updated the time if not early stopping
      {if (l<no_cohort)
        time<- time+max(V_T,V_E)
      else if (l==no_cohort)
        time=time+min(max(temp_D[,2],temp_D[,5]),max(U_T,U_E))
      }
    }
    
    dose_history[l+1]<-dose
    
  }
  ##### Summary the result
  ###### Final dose selection
  Big_summary[k,1]<-dose_history[length(dose_history)]
  
  ##### Get early stopping
  
  ###### Patient allocation gor dose 1-5
  Big_summary[k,2]<-sum(accrue_D[,1]==1)
  Big_summary[k,3]<-sum(accrue_D[,1]==2)
  Big_summary[k,4]<-sum(accrue_D[,1]==3)
  Big_summary[k,5]<-sum(accrue_D[,1]==4)
  Big_summary[k,6]<-sum(accrue_D[,1]==5)
  
  ##### Plus total toxiciy subject number and total non-efficacy (futility) event number 
  Big_summary[k,7]<-sum(accrue_D[,4])
  Big_summary[k,8]<-sum(accrue_D[,7])
  
  ##### Get early stopping
  Big_summary[k,9]<-early_stopping
  
  ##### Get the trial time
  Big_summary[k,10]<- time
  
  ##### Get the stage I number
  Big_summary[k,11]<- sum(stage==1)
  
  ##### Get the stage II number
  Big_summary[k,12]<- sum(stage==2)
  
  ##### Get the dose history
  dose_history_sum[k,]<- dose_history
}

#########
dose_finding<-numeric(length = nrow(para_U)+1)
 
 for (i in 1:length(dose_finding)) { 
dose_finding[i]<-sum(Big_summary[,1]==i-(1))
}

summ<-apply(Big_summary[,2:12],2,mean)

return(list(dose_finding,summ))

}


###### U-BOIN Eameple 

###### scenario 5 input parameters

DLT<- c(0.1,0.3,0.5,0.55,0.65)

efficacy<- c(0.45,0.45,0.45,0.45,0.45)

###### Output weibull parameters for further simulation
para_BOIN<-para_finding_BOIN(U_T=3,V_T=2,U_E=3,V_E=2,event_rate=0.6, DLT=DLT,eff=efficacy)

################
set.seed(333333)

result_BOIN<-CWL_U_BOIN_OC(para_U=para_BOIN,phi_0=0.3,U_T=3,U_E=3,V_T=2,V_E=2,p_cut=0.95,q_cut=0.9,
                   pi_T=0.3,pi_E=0.2, no_cohort=20, cohort_size=3, score=c(0,30,50,100),M=1000)

##############################################################################################
##############################################################################################
########### The following part is for CWL-OR function

######## We use following function to get the suitable parameter of weibull distribution for OR method simualtion
#' @param U_T is the assessment period for DLT
#' @param v_T is the interim period for DLT
#' @param U_E is the assessment period for efficacy
#' @param v_E is the interim period for efficacy
#' @param DLT is DLT rate curve for each dose level within assessment period
#' @param eff is efficacy rate curve for each dose level within assessment period
#' @param efficacy_0 is the efficacy rate curve condition on non-DLT outcome
#' @param event_rate is the efficacy and DLT event rate happens within V_T & V_E

para_find_OR<-function(U_T=3,V_T=2,U_E=3,V_E=2, event_rate=0.6, DLT=DLT,eff=efficacy, eff_0= efficacy_0) {
  
  para_U<- matrix(0,nrow=5,ncol=6)
  colnames(para_U)<-c("p_k","p_lambda","q0_k","q0_lambda","q1_k","q1_lambda")
  
  for (i in 1:length(DLT)) {
    
    ###### DLT rate
    p_U<-DLT[i]
    p_Z<-DLT[i]*event_rate
    
    
    ###### Efficacy_condition non_DLT
    q0_U<-efficacy_0[i]
    #######
    q0_Z<-q0_U*event_rate
    
    ###### Efficacy_condition DLT
    q1_U<-(efficacy[i]-efficacy_0[i]*(1-DLT[i]))/DLT[i]
    #######
    q1_Z<-q1_U*event_rate
    
    #### pars[1]=shape=k pars[2]=lambda=scale
    f_1<-function(pars){
      
      y<- numeric(2)
      y[1]<-pweibull(U_T, pars[1],pars[2], lower.tail = TRUE, log.p = FALSE) - p_U
      y[2]<-pweibull(V_T, pars[1],pars[2], lower.tail = TRUE, log.p = FALSE) - p_Z
      
      return (y)
    }
    
    
    ### xstart is start point , nleqslv is the function to find parameter
    xstart <- c(2,3)
    
    outcome1<-nleqslv(xstart, fn=f_1, method="Newton",global="qline",control = list(allowSingular=TRUE))
    
    
    #### pars[1]=shape=k_q0 pars[2]=scale=lambda_q0
    f_2<-function(pars){
      
      y<- numeric(2)
      y[1]<-pweibull(U_E, pars[1],pars[2], lower.tail = TRUE, log.p = FALSE) - q0_U
      y[2]<-pweibull(V_E, pars[1],pars[2], lower.tail = TRUE, log.p = FALSE) - q0_Z
      
      return (y)
    }
    
    #### pars[1]=shape=k_q1 pars[2]=scale=lambda_q1
    f_3<-function(pars){
      
      y<- numeric(2)
      y[1]<-pweibull(U_E, pars[1],pars[2], lower.tail = TRUE, log.p = FALSE) - q1_U
      y[2]<-pweibull(V_E, pars[1],pars[2], lower.tail = TRUE, log.p = FALSE) - q1_Z
      
      return (y)
    }
    
    xstart<-c(2,3)
    
    outcome2<-nleqslv(xstart, fn=f_2, method="Newton",global="qline",control = list(allowSingular=TRUE))
    outcome3<-nleqslv(xstart, fn=f_3, method="Newton",global="qline",control = list(allowSingular=TRUE))
    
    
    ### get the result for toxiciy's parameters k and lambda
    k<-outcome1$x[1]
    lambda<-outcome1$x[2]
    
    
    ### get the result for conditional efficacy's parameters k and lambda
    k_q0<-outcome2$x[1]
    lambda_q0<-outcome2$x[2]
    
    k_q1<-outcome3$x[1]
    lambda_q1<-outcome3$x[2]
    
    ###### Save the para
    para_U[i,]<-c(k,lambda,k_q0,lambda_q0,k_q1,lambda_q1)
    
  }
  
  return(para_U)
  
}


######### We use folloing function to do the simulation for CWL OR method
#' @param para is the initial parameter result geeting from above function 
#' @param U_T is the assessment period for DLT
#' @param v_T is the interim period for DLT
#' @param U_E is the assessment period for efficacy
#' @param v_E is the interim period for efficacy
#' @param pi_T is pre-setting DLT rate upper boundary for the medicine
#' @param pi_E is pre-setting efficacy rate lower boundary for the medicine
#' @param p_cut is the posterior distribution cut-off for DLT rate
#' @param p_escl is the posterior distribution DLT cut-off for dose climbing
#' @param q_cut is the posterior distribution cut-off for efficacy rate
#' @param no_cohort is the totol cohorts involved in one clinicl trial
#' @param cohort_size is the cohort size for one cohort
#' @param M is the total simualtion times.

CWL_OR_OC<-function(para_U=para,U_T=3,U_E=3,V_T=2,V_E=2,p_cut=0.9,p_escl=0.3,q_cut=0.9, 
                    pi_T=0.3,pi_E=0.25, no_cohort=20,cohort_size=3, M=1000) {
  
  #### Write the model with 5 dose levels
  model_5d<- "model{

#   # Likelihood
y1~dmulti(R1,n1)
y2~dmulti(R2,n2)
y3~dmulti(R3,n3)
y4~dmulti(R4,n4)
y5~dmulti(R5,n5)
y6~dmulti(R6,n6)
y7~dmulti(R7,n7)
y8~dmulti(R8,n8)
y9~dmulti(R9,n9)
y10~dmulti(R10,n10)



#   # Prior
#   # 1 is toxicity and efficacy; 2 is toxicity and non-efficacy
#   # 3 is non-toxiciy and efficacy;4 is non-toxicity and non-efficacy;
##### Fully observed
R1[1]<- p_1*q_11
R1[2]<- p_1*(1-q_11)
R1[3]<- (1-p_1)*q_01
R1[4]<- (1-p_1)*(1-q_01)


R2[1]<- p_2*q_12
R2[2]<- p_2*(1-q_12)
R2[3]<- (1-p_2)*q_02
R2[4]<- (1-p_2)*(1-q_02)


R3[1]<- p_3*q_13
R3[2]<- p_3*(1-q_13)
R3[3]<- (1-p_3)*q_03
R3[4]<- (1-p_3)*(1-q_03)


R4[1]<- p_4*q_14
R4[2]<- p_4*(1-q_14)
R4[3]<- (1-p_4)*q_04
R4[4]<- (1-p_4)*(1-q_04)


R5[1]<- p_5*q_15
R5[2]<- p_5*(1-q_15)
R5[3]<- (1-p_5)*q_05
R5[4]<- (1-p_5)*(1-q_05)


##### Partly missing----2 tau
R6[1]<- (V_E/U_E)*q_11*(V_T/U_T)*p_1
R6[2]<- (V_T/U_T)*p_1*(1-(V_E/U_E)*q_11)
R6[3]<- (1-(V_T/U_T)*p_1)*(V_E/U_E)*(q_11*((1-V_T/U_T)*p_1/((1-V_T/U_T)*p_1+1-p_1))+q_01*((1-p_1)/((1-V_T/U_T)*p_1+1-p_1)))
R6[4]<- (1-(V_T/U_T)*p_1)*(1-(V_E/U_E)*(q_11*((1-V_T/U_T)*p_1/((1-V_T/U_T)*p_1+1-p_1))+q_01*((1-p_1)/((1-V_T/U_T)*p_1+1-p_1))))


R7[1]<- (V_E/U_E)*q_12*(V_T/U_T)*p_2
R7[2]<- (V_T/U_T)*p_2*(1-(V_E/U_E)*q_12)
R7[3]<- (1-(V_T/U_T)*p_2)*(V_E/U_E)*(q_12*((1-V_T/U_T)*p_2/((1-V_T/U_T)*p_2+1-p_2))+q_02*((1-p_2)/((1-V_T/U_T)*p_2+1-p_2)))
R7[4]<- (1-(V_T/U_T)*p_2)*(1-(V_E/U_E)*(q_12*((1-V_T/U_T)*p_2/((1-V_T/U_T)*p_2+1-p_2))+q_02*((1-p_2)/((1-V_T/U_T)*p_2+1-p_2))))


R8[1]<- (V_E/U_E)*q_13*(V_T/U_T)*p_3
R8[2]<- (V_T/U_T)*p_3*(1-(V_E/U_E)*q_13)
R8[3]<- (1-(V_T/U_T)*p_3)*(V_E/U_E)*(q_13*((1-V_T/U_T)*p_3/((1-V_T/U_T)*p_3+1-p_3))+q_03*((1-p_3)/((1-V_T/U_T)*p_3+1-p_3)))
R8[4]<- (1-(V_T/U_T)*p_3)*(1-(V_E/U_E)*(q_13*((1-V_T/U_T)*p_3/((1-V_T/U_T)*p_3+1-p_3))+q_03*((1-p_3)/((1-V_T/U_T)*p_3+1-p_3))))


R9[1]<- (V_E/U_E)*q_14*(V_T/U_T)*p_4
R9[2]<- (V_T/U_T)*p_4*(1-(V_E/U_E)*q_14)
R9[3]<- (1-(V_T/U_T)*p_4)*(V_E/U_E)*(q_14*((1-V_T/U_T)*p_4/((1-V_T/U_T)*p_4+1-p_4))+q_04*((1-p_4)/((1-V_T/U_T)*p_4+1-p_4)))
R9[4]<- (1-(V_T/U_T)*p_4)*(1-(V_E/U_E)*(q_14*((1-V_T/U_T)*p_4/((1-V_T/U_T)*p_4+1-p_4))+q_04*((1-p_4)/((1-V_T/U_T)*p_4+1-p_4))))


R10[1]<- (V_E/U_E)*q_15*(V_T/U_T)*p_5
R10[2]<- (V_T/U_T)*p_5*(1-(V_E/U_E)*q_15)
R10[3]<- (1-(V_T/U_T)*p_5)*(V_E/U_E)*(q_15*((1-V_T/U_T)*p_5/((1-V_T/U_T)*p_5+1-p_5))+q_05*((1-p_5)/((1-V_T/U_T)*p_5+1-p_5)))
R10[4]<- (1-(V_T/U_T)*p_5)*(1-(V_E/U_E)*(q_15*((1-V_T/U_T)*p_5/((1-V_T/U_T)*p_5+1-p_5))+q_05*((1-p_5)/((1-V_T/U_T)*p_5+1-p_5))))

######## Toxicity and efficacy planned follow-up time and fixed follow-up time

U_E<-3
U_T<-3

####### actual follow-up time
V_E<-2
V_T<-2


########
### Dose 1
p_1<- exp(lambda_1)/(1+exp(lambda_1))
q_01<-exp(mu_1)/(1+exp(mu_1))
q_11<-exp(nu_1)/(1+exp(nu_1))


### Dose 2
p_2<- (exp(lambda_1)+ exp(lambda_2))/(1+exp(lambda_1)+ exp(lambda_2))
q_02<-exp(mu_1+mu_2)/(1+exp(mu_1+mu_2))
q_12<-exp(nu_1+nu_2)/(1+exp(nu_1+nu_2))


### Dose 3
p_3<- (exp(lambda_1)+ exp(lambda_2)+ exp(lambda_3))/(1+ exp(lambda_1)+ exp(lambda_2)+ exp(lambda_3))
q_03<-exp(mu_1+ mu_2+ mu_3)/(1+exp(mu_1+ mu_2+ mu_3))
q_13<-exp(nu_1+ nu_2+ nu_3)/(1+exp(nu_1+ nu_2+ nu_3))

### Dose 4
p_4<- (exp(lambda_1)+ exp(lambda_2)+ exp(lambda_3)+ exp(lambda_4))/(1+ exp(lambda_1)+ exp(lambda_2)+ exp(lambda_3)+ exp(lambda_4))
q_04<-exp(mu_1+ mu_2+ mu_3+ mu_4)/(1+exp(mu_1+ mu_2+ mu_3+ mu_4))
q_14<-exp(nu_1+ nu_2+ nu_3+ nu_4)/(1+exp(nu_1+ nu_2+ nu_3+ nu_4))

### Dose 5
p_5<- (exp(lambda_1)+ exp(lambda_2)+ exp(lambda_3)+ exp(lambda_4)+ exp(lambda_5))/(1+ exp(lambda_1)+ exp(lambda_2)+ exp(lambda_3)+ exp(lambda_4)+ exp(lambda_5))
q_05<-exp(mu_1+ mu_2+ mu_3+ mu_4+ mu_5)/(1+exp(mu_1+ mu_2+ mu_3+ mu_4+ mu_5))
q_15<-exp(nu_1+ nu_2+ nu_3+ nu_4+ nu_5)/(1+exp(nu_1+ nu_2+ nu_3+ nu_4+ nu_5))

####### Define prior for all 15 parameters
lambda_1~dnorm(-3,0.01)
lambda_2~dnorm(-3,0.01)
lambda_3~dnorm(-3,0.01)
lambda_4~dnorm(-3,0.01)
lambda_5~dnorm(-3,0.01)

mu_1~dnorm(0,0.01)
mu_2~dnorm(0,0.01)
mu_3~dnorm(0,0.01)
mu_4~dnorm(0,0.01)
mu_5~dnorm(0,0.01)

nu_1~dnorm(0,0.01)
nu_2~dnorm(0,0.01)
nu_3~dnorm(0,0.01)
nu_4~dnorm(0,0.01)
nu_5~dnorm(0,0.01)

}"
  
  
  ###################  
  options(digits=5)
  
  #### Create the Big summary table
  Big_summary<-matrix(0,M,10)
  
  dose_history_sum<-matrix(0,M,(no_cohort+1))
  
  for (k in 1:M) {
    ### create the accrue data frame, accrue_s_a, accrue_delta_a
    accrue_D<-matrix(8,1,7)
    accrue_D<-accrue_D[-1,]
    
    #### Rest the parameter
    dose<-1
    
    time<-0
    
    early_stopping <- 0
    
    sum_sim_result<-matrix(0,5,4)
    ##### The last one is final choice dose, not given to subjects
    dose_history<-c(1,rep(0,no_cohort))
    
    for (l in 1:no_cohort){
      
      #### Tempory data
      temp_D<-matrix(0,cohort_size,7)
      
      temp_D[,1]<- rep(dose,cohort_size)
      
      for (i in 1:cohort_size){
        
        temp_D[i,2]<-rweibull(1,shape = para_U[dose,1], scale = para_U[dose,2])
        temp_D[i,3]<-(temp_D[i,2]<=V_T)
        temp_D[i,4]<-(temp_D[i,2]<=U_T)
        
        #### Based on toxicity to simulate efficacy  
        temp_D[i,5]<-ifelse((temp_D[i,2]>=U_T),rweibull(1,shape = para_U[dose,3], scale = para_U[dose,4]),
                            rweibull(1,shape = para_U[dose,5], scale = para_U[dose,6]))
        
        temp_D[i,6]<-(temp_D[i,5]<=V_E)
        temp_D[i,7]<-(temp_D[i,5]<=U_E)
        
      }
      
      accrue_D<-rbind(accrue_D,temp_D)
      
      #### Compute y according to accrue_D
      #### y1-y5 Fully observed ,y6-y10 is partly missing
      y<-matrix(0,10,4)
      
      if (l==1)
      {y[1,1]<-sum(accrue_D[,4]*accrue_D[,7])
      y[1,2]<-sum(accrue_D[,4]*(1-accrue_D[,7]))
      y[1,3]<-sum((1-accrue_D[,4])*accrue_D[,7])
      y[1,4]<-sum((1-accrue_D[,4])*(1-accrue_D[,7]))
      } 
      else if (l<no_cohort)
      { for (j in 1:5){
        y[j,1]<-sum(accrue_D[1:(cohort_size*l-cohort_size),4]*accrue_D[1:(cohort_size*l-cohort_size),7]*(accrue_D[1:(cohort_size*l-cohort_size),1]==j))
        y[j,2]<-sum(accrue_D[1:(cohort_size*l-cohort_size),4]*(1-accrue_D[1:(cohort_size*l-cohort_size),7])*(accrue_D[1:(cohort_size*l-cohort_size),1]==j))
        y[j,3]<-sum((1-accrue_D[1:(cohort_size*l-cohort_size),4])*accrue_D[1:(cohort_size*l-cohort_size),7]*(accrue_D[1:(cohort_size*l-cohort_size),1]==j))
        y[j,4]<-sum((1-accrue_D[1:(cohort_size*l-cohort_size),4])*(1-accrue_D[1:(cohort_size*l-cohort_size),7])*(accrue_D[1:(cohort_size*l-cohort_size),1]==j))
      }
        for (j in 6:10){
          y[j,1]<-sum(accrue_D[(cohort_size*l-2):(cohort_size*l),3]*accrue_D[(cohort_size*l-2):(cohort_size*l),6]*(accrue_D[(cohort_size*l-2):(cohort_size*l),1]==(j-5)))
          y[j,2]<-sum(accrue_D[(cohort_size*l-2):(cohort_size*l),3]*(1-accrue_D[(cohort_size*l-2):(cohort_size*l),6])*(accrue_D[(cohort_size*l-2):(cohort_size*l),1]==(j-5)))
          y[j,3]<-sum((1-accrue_D[(cohort_size*l-2):(cohort_size*l),3])*accrue_D[(cohort_size*l-2):(cohort_size*l),6]*(accrue_D[(cohort_size*l-2):(cohort_size*l),1]==(j-5)))
          y[j,4]<-sum((1-accrue_D[(cohort_size*l-2):(cohort_size*l),3])*(1-accrue_D[(cohort_size*l-2):(cohort_size*l),6])*(accrue_D[(cohort_size*l-2):(cohort_size*l),1]==(j-5)))
        }
      }
      ######
      else if (l==no_cohort)
      {for (j in 1:5){
        y[j,1]<-sum(accrue_D[,4]*accrue_D[,7]*(accrue_D[,1]==j))
        y[j,2]<-sum(accrue_D[,4]*(1-accrue_D[,7])*(accrue_D[,1]==j))
        y[j,3]<-sum((1-accrue_D[,4])*accrue_D[,7]*(accrue_D[,1]==j))
        y[j,4]<-sum((1-accrue_D[,4])*(1-accrue_D[,7])*(accrue_D[,1]==j))
      }
      }
      
      y1<-y[1,]
      y2<-y[2,]
      y3<-y[3,]
      y4<-y[4,]
      y5<-y[5,]
      y6<-y[6,]
      y7<-y[7,]
      y8<-y[8,]
      y9<-y[9,]
      y10<-y[10,]
      #####
      n1<-sum(y1)
      n2<-sum(y2)
      n3<-sum(y3)
      n4<-sum(y4)
      n5<-sum(y5)
      #####
      n6<-sum(y6)
      n7<-sum(y7)
      n8<-sum(y8)
      n9<-sum(y9)
      n10<-sum(y10)
      
      ######
      model <- jags.model(textConnection(model_5d), data = list(y1=y1,y2=y2,y3=y3,y4=y4,y5=y5,y6=y6,y7=y7,y8=y8,y9=y9,y10=y10, n1=n1,n2=n2,n3=n3,n4=n4,n5=n5,n6=n6,n7=n7,n8=n8,n9=n9,n10=n10),quiet=T)
      
      update(model, 1000, progress.bar="none")
      
      S<-5000
      
      samp <- jags.samples(model, n.iter=S, variable.names=c("lambda_1","lambda_2","lambda_3","lambda_4","lambda_5","mu_1","mu_2","mu_3","mu_4","mu_5","nu_1","nu_2","nu_3","nu_4","nu_5"), 
                           progress.bar="none")
      
      samp.lambda_1<- samp$lambda_1[1:S]
      samp.lambda_2<- samp$lambda_2[1:S]
      samp.lambda_3<- samp$lambda_3[1:S]
      samp.lambda_4<- samp$lambda_4[1:S]
      samp.lambda_5<- samp$lambda_5[1:S]
      
      samp.mu_1<- samp$mu_1[1:S]
      samp.mu_2<- samp$mu_2[1:S]
      samp.mu_3<- samp$mu_3[1:S]
      samp.mu_4<- samp$mu_4[1:S]
      samp.mu_5<- samp$mu_5[1:S]
      
      samp.nu_1<- samp$nu_1[1:S]
      samp.nu_2<- samp$nu_2[1:S]
      samp.nu_3<- samp$nu_3[1:S]
      samp.nu_4<- samp$nu_4[1:S]
      samp.nu_5<- samp$nu_5[1:S]
      
      
      #### Get the result from p1-p5, q01-q05, q11-q15;
      #### Dose 1
      p_1<-exp(samp.lambda_1)/(1+exp(samp.lambda_1))
      q_01<-exp(samp.mu_1)/(1+exp(samp.mu_1))
      q_11<-exp(samp.nu_1)/(1+exp(samp.nu_1))
      q_1<- p_1*q_11+(1-p_1)*q_01
      
      ### pi_1_01<- (1-p_1)*q_01
      ### pi_1_00<- (1-p_1)-pi_1_01
      
      ### Dose 2
      p_2<- (exp(samp.lambda_1)+ exp(samp.lambda_2))/(1+exp(samp.lambda_1)+ exp(samp.lambda_2))
      q_02<-exp(samp.mu_1+samp.mu_2)/(1+exp(samp.mu_1+samp.mu_2))
      q_12<-exp(samp.nu_1+samp.nu_2)/(1+exp(samp.nu_1+samp.nu_2))
      q_2<- p_2*q_12+(1-p_2)*q_02
      
      ### pi_2_01<- (1-p_2)*q_02
      ### pi_2_00<- (1-p_2)-pi_2_01
      
      ### Dose 3
      p_3<- (exp(samp.lambda_1)+ exp(samp.lambda_2)+ exp(samp.lambda_3))/(1+ exp(samp.lambda_1)+ exp(samp.lambda_2)+ exp(samp.lambda_3))
      q_03<-exp(samp.mu_1+ samp.mu_2+ samp.mu_3)/(1+exp(samp.mu_1+ samp.mu_2+ samp.mu_3))
      q_13<-exp(samp.nu_1+ samp.nu_2+ samp.nu_3)/(1+exp(samp.nu_1+ samp.nu_2+ samp.nu_3))
      q_3<- p_3*q_13+(1-p_3)*q_03
      
      ### pi_3_01<- (1-p_3)*q_03
      ### pi_3_00<- (1-p_3)-pi_3_01
      
      ### Dose 4
      p_4<- (exp(samp.lambda_1)+ exp(samp.lambda_2)+ exp(samp.lambda_3)+ exp(samp.lambda_4))/(1+ 
                                                                                                exp(samp.lambda_1)+ exp(samp.lambda_2)+ exp(samp.lambda_3)+ exp(samp.lambda_4))
      q_04<-exp(samp.mu_1+ samp.mu_2+ samp.mu_3+ samp.mu_4)/(1+
                                                               exp(samp.mu_1+ samp.mu_2+ samp.mu_3+ samp.mu_4))
      q_14<-exp(samp.nu_1+ samp.nu_2+ samp.nu_3+ samp.nu_4)/(1+
                                                               exp(samp.nu_1+ samp.nu_2+ samp.nu_3+ samp.nu_4))
      q_4<- p_4*q_14+(1-p_4)*q_04
      
      ### pi_4_01<- (1-p_4)*q_04
      ### pi_4_00<- (1-p_4)-pi_4_01
      
      
      ### Dose 5
      p_5<- (exp(samp.lambda_1)+ exp(samp.lambda_2)+ exp(samp.lambda_3)+ exp(samp.lambda_4)+ exp(samp.lambda_5))/(1+
                                                                                                                    exp(samp.lambda_1)+ exp(samp.lambda_2)+ exp(samp.lambda_3)+ exp(samp.lambda_4)+ exp(samp.lambda_5))
      q_05<-exp(samp.mu_1+ samp.mu_2+ samp.mu_3+ samp.mu_4+ samp.mu_5)/(1+
                                                                          exp(samp.mu_1+ samp.mu_2+ samp.mu_3+ samp.mu_4+ samp.mu_5))
      q_15<-exp(samp.nu_1+ samp.nu_2+ samp.nu_3+ samp.nu_4+ samp.nu_5)/(1+
                                                                          exp(samp.nu_1+ samp.nu_2+ samp.nu_3+ samp.nu_4+ samp.nu_5))
      q_5<- p_5*q_15+(1-p_5)*q_05
      
      ### pi_5_01<- (1-p_5)*q_05
      ### pi_5_00<- (1-p_5)-pi_5_01
      
      ##### Combine the original result
      sim_result<-cbind(p_1,q_1,q_01,q_11,p_2,q_2,q_02,q_12,
                        p_3,q_3,q_03,q_13,p_4,q_4,q_04,q_14,p_5,q_5,q_05,q_15)
      
      
      for (j in 1:5) {
        
        sum_sim_result[j,]<- apply(sim_result[ ,((j-1)*4+1):((j-1)*4+4)],2,mean)
        
      }
      
      
      ##### Get accpetable dose level set
      acceptable_dose<-rep(0,5)
      
      for (j in 1:5){
        
        acceptable_dose[j]<-(((sum(sim_result[,1+(j-1)*4]<pi_T)/S)>(1-p_cut))*((sum(sim_result[,2+(j-1)*4]>pi_E)/S)>(1-q_cut)))
        
      }
      
      ##### 
      omega_2<-omega_3<-rep(0,5)
      
      for (j in 1:5){
        
        omega_2[j]<-mean(sim_result[,1+(j-1)*4])*(1-mean(sim_result[,2+(j-1)*4]))/((1-mean(sim_result[,1+(j-1)*4]))*mean(sim_result[,2+(j-1)*4]))
        
        omega_3[j]<- omega_2[j]*(1-mean(sim_result[,3+(j-1)*4]))/mean(sim_result[,3+(j-1)*4])
      }
      
      ######
      ##### This code without stopping rule and time account
      
      max_dose<-max(dose_history)
      
      {if ((max_dose<=4)*((sum(sim_result[,(max_dose-1)*4+1]<pi_T)/S)>p_escl))
        
        dose<-max_dose+1
        
        else if (sum(acceptable_dose)==0)
        {dose<-0
        
        if (l==1)  {
          early_stopping <- 1
          time<-max(U_T,U_E)
        } 
        else if  (l==no_cohort) {
          early_stopping <- 0
          time<- time+min(max(temp_D[,2],temp_D[,5]),max(U_T,U_E))
        }
        else if  (l<no_cohort) {
          early_stopping <- 1
          time<- time+min(max(temp_D[,2],temp_D[,5]),max(V_T,V_E))
        }
        break
        }
        
        
        else if (sum(acceptable_dose)>0)
          
        {
          min_value<-min(subset(omega_3,acceptable_dose*omega_3!=0))
          
          dose<-which(omega_3==min_value)
        }
        
      }
      
      ##### Updated the time if not early stopping
      {if (l==1)
        time<-max(U_T,U_E)
        
        else if (l==no_cohort)
          
          time=time+min(max(temp_D[,2],temp_D[,5]),max(U_T,U_E))
        
        else if ((l>1)*(l<no_cohort))
          
          time=time+max(V_T,V_E)
      }
      
      dose_history[l+1]<-dose
      
    }
    
    ###### Final dose selection
    Big_summary[k,1]<-dose_history[21]
    
    ###### Patient allocation gor dose 1-5
    Big_summary[k,2]<-sum(accrue_D[,1]==1)
    Big_summary[k,3]<-sum(accrue_D[,1]==2)
    Big_summary[k,4]<-sum(accrue_D[,1]==3)
    Big_summary[k,5]<-sum(accrue_D[,1]==4)
    Big_summary[k,6]<-sum(accrue_D[,1]==5)
    
    ##### Plus total toxiciy subject number and total non-efficacy (futility) event number 
    Big_summary[k,7]<-sum(accrue_D[,4])
    Big_summary[k,8]<-(nrow(accrue_D)-sum(accrue_D[,7]))
    
    ##### Get early stopping
    Big_summary[k,9]<-early_stopping
    
    ##### Get the trial time
    Big_summary[k,10]<- time
    
    ##### Get the dose history
    dose_history_sum[k, ]<- dose_history
    
  }  
  
  #########
  dose_finding<-numeric(length = nrow(para_U)+1)
  
  for (i in 1:length(dose_finding)) { 
    dose_finding[i]<-sum(Big_summary[,1]==i-(1))
  }
  
  summ<-apply(Big_summary[,2:10],2,mean)
  
  return(list(dose_finding,summ))
  
}


#### Example OR method 
##### Scenario 1

DLT<-c(0.05,0.25,0.28,0.35,0.5)
##### efficacy
efficacy<-c(0.4,0.45,0.45,0.4,0.3)
##### efficacy condition on Non-DLT
efficacy_0<-c(0.4,0.3,0.25,0.25,0.25)

######  Output weibull parameters for further simulation
para_OR<-para_find_OR(U_T=3,V_T=2,U_E=3,V_E=2, event_rate=0.6, DLT=DLT,eff=efficacy, eff_0= efficacy_0)

############ DO the Simulation for OR method

resut_OR<-CWL_OR_OC(para_U=para_OR,U_T=3,U_E=3,V_T=2,V_E=2,p_cut=0.9,p_escl=0.3,q_cut=0.9, 
                    pi_T=0.3,pi_E=0.25, no_cohort=20,cohort_size=3, M=1000)



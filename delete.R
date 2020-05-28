Vn = c("Time0", "MMFdos0", "age", "Sex", "ACA", "SCL70", "RNAPol", "cfcbFVC0q", "cfcbmRSS0q", "Race_1", "Race_2", "ethnic_0", "ethnic_1", "YTime")

V =  c("Time0", "MMFdos0", "age", "Sex", "ACA", "SCL70", "RNAPol", "cfcbFVC0q", "cfcbmRSS0q", "Race_1", "Race_2", "ethnic_0", "ethnic_1")

(FVC1q, mRSS1q, A1) ~ (FVC0, mRSS0, V)

V = MMFdos0 + age + Sex + ACA + SCL70 + RNAPol + Race_1 + Race_2 + ethnic_0 + ethnic_1
Z = Time0 + YTime

  
  Y1q ~  A1 + V (+ Zb)
  
  A1  ~  Y0 + A0 + V (+ b)
  
  Y0 = FVC0 + mRSS0
       FVC0q + mRSS0q
       cfcbFVC0 + cfcbmRSS0
       cfcbFVC0q + cfcbmRSS0q
       
  MCMCglmm(
    fixed = cbind(FVC1q, mRSS1q, A1) , data = dat)
  
  
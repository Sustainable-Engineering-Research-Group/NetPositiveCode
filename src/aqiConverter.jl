function InvLinear(AQIhigh, AQIlow, Conchigh, Conclow, a)
    c=((a-AQIlow)/(AQIhigh-AQIlow))*(Conchigh-Conclow)+Conclow
    return c
end

function ConcPM25(a)
    ismissing(a) && return a
    if (a>=0 && a<=50)
        ConcCalc=InvLinear(50,0,12,0,a)
    elseif (a>50 && a≤100)
        ConcCalc=InvLinear(100,51,35.4,12.1,a)
    elseif (a>100 && a≤150)
        ConcCalc=InvLinear(150,101,55.4,35.5,a)
    elseif (a>150 && a≤200)
        ConcCalc=InvLinear(200,151,150.4,55.5,a)
    elseif (a>200 && a<=300)
        ConcCalc=InvLinear(300,201,250.4,150.5,a)
    elseif (a>300 && a<=400)
        ConcCalc=InvLinear(400,301,350.4,250.5,a)
    elseif (a>400 && a<=500)
        ConcCalc=InvLinear(500,401,500.4,350.5,a)
    else
        ConcCalc=missing
    end
    return ConcCalc
end

function ConcPM10(a)
    ismissing(a) && return a
    if (a>=0 && a<=50)
        ConcCalc=InvLinear(50,0,54,0,a)
    elseif (a>50 && a<=100)
        ConcCalc=InvLinear(100,51,154,55,a)
    elseif (a>100 && a<=150)
        ConcCalc=InvLinear(150,101,254,155,a)
    elseif (a>150 && a<=200)
        ConcCalc=InvLinear(200,151,354,255,a)
    elseif (a>200 && a<=300)
        ConcCalc=InvLinear(300,201,424,355,a)
    elseif (a>300 && a<=400)
        ConcCalc=InvLinear(400,301,504,425,a)
    elseif (a>400 && a<=500)
        ConcCalc=InvLinear(500,401,604,505,a)
    else
        ConcCalc=missing
    end
    return ConcCalc
end

function ConcCO(a)
    ismissing(a) && return a
    if (a>=0 && a<=50)
        ConcCalc=InvLinear(50,0,4.4,0,a)
    elseif (a>50 && a<=100)
        ConcCalc=InvLinear(100,51,9.4,4.5,a)
    elseif (a>100 && a<=150)
        ConcCalc=InvLinear(150,101,12.4,9.5,a)
    elseif (a>150 && a<=200)
        ConcCalc=InvLinear(200,151,15.4,12.5,a)
    elseif (a>200 && a<=300)
        ConcCalc=InvLinear(300,201,30.4,15.5,a)
    elseif (a>300 && a<=400)
        ConcCalc=InvLinear(400,301,40.4,30.5,a)
    elseif (a>400 && a<=500)
        ConcCalc=InvLinear(500,401,50.4,40.5,a)
    else
        ConcCalc=missing
    end
    return ConcCalc/0.873e-3
end

function ConcSO21hr(a)
    ismissing(a) && return a
    if (a>=0 && a<=50)
        ConcCalc=InvLinear(50,0,35,0,a)
    elseif (a>50 && a<=100)
        ConcCalc=InvLinear(100,51,75,36,a)
    elseif (a>100 && a<=150)
        ConcCalc=InvLinear(150,101,185,76,a)
    elseif (a>150 && a<=200)
        ConcCalc=InvLinear(200,151,304,186,a)
    else
        ConcCalc=missing
    end
    return ConcCalc/0.382
end

function ConcSO224hr(a)
    ismissing(a) && return a
    if (a<=200)
        ConcCalc=304
    elseif (a>=201 && a<=300)
        ConcCalc=InvLinear(300,201,604,305,a)
    elseif (a>300 && a<=400)
        ConcCalc=InvLinear(400,301,804,605,a)
    elseif (a>400 && a<=500)
        ConcCalc=InvLinear(500,401,1004,805,a)
    else
        ConcCalc=missing
    end
    return ConcCalc/0.382
end

function ConcOzone8hr(a)
    ismissing(a) && return a
    if (a>=0 && a<=50)
        ConcCalc=InvLinear(50,0,54,0,a)
    elseif (a>50 && a<=100)
        ConcCalc=InvLinear(100,51,70,55,a)
    elseif (a>100 && a<=150)
        ConcCalc=InvLinear(150,101,85,71,a)
    elseif (a>150 && a<=200)
        ConcCalc=InvLinear(200,151,105,86,a)
    elseif (a>200 && a<=300)
        ConcCalc=InvLinear(300,201,200,106,a)
    else
        ConcCalc=missing
    end
    return ConcCalc/0.509
end

function ConcOzone1hr(a)
    ismissing(a) && return a
    if a≤100
        ConcCalc = 124
    elseif (a>100 && a<=150)
        ConcCalc=InvLinear(150,101,164,125,a)
    elseif (a>150 && a<=200)
        ConcCalc=InvLinear(200,151,204,165,a)
    elseif (a>200 && a<=300)
        ConcCalc=InvLinear(300,201,404,205,a)
    elseif (a>300 && a<=400)
        ConcCalc=InvLinear(400,301,504,405,a)
    elseif (a>400 && a<=500)
        ConcCalc=InvLinear(500,401,604,505,a)
    else
        ConcCalc=missing
    end
    return ConcCalc/0.509
end

function ConcNO2(a)
    ismissing(a) && return a
    if (a>=0 && a<=50)
        ConcCalc=InvLinear(50,0,53,0,a)
    elseif (a>50 && a<=100)
        ConcCalc=InvLinear(100,51,100,54,a)
    elseif (a>100 && a<=150)
        ConcCalc=InvLinear(150,101,360,101,a)
    elseif (a>150 && a<=200)
        ConcCalc=InvLinear(200,151,649,361,a)
    elseif (a>200 && a<=300)
        ConcCalc=InvLinear(300,201,1244,650,a)
    elseif (a>300 && a<=400)
        ConcCalc=InvLinear(400,301,1644,1245,a)
    elseif (a>400 && a<=500)
        ConcCalc=InvLinear(500,401,2044,1645,a)
    else
        ConcCalc=missing
    end
    return ConcCalc/0.532
end

function Linear(AQIhigh, AQIlow, Conchigh, Conclow, Concentration)
    Conc = Concentration
    a = ((Conc-Conclow)/(Conchigh-Conclow))*(AQIhigh-AQIlow)+AQIlow
    linear = ismissing(a) ? missing : round(Int, a)
    return linear
end

function AQIPM25(Concentration)
    Conc = Concentration
    ismissing(Conc) && return missing
    c = round(Conc; digits = 1)
    if (c>=0 && c<12.1)
        AQI=Linear(50,0,12,0,c)
    elseif (c>=12.1 && c<35.5)
        AQI=Linear(100,51,35.4,12.1,c)
    elseif (c>=35.5 && c<55.5)
        AQI=Linear(150,101,55.4,35.5,c)
    elseif (c>=55.5 && c<150.5)
        AQI=Linear(200,151,150.4,55.5,c)
    elseif (c>=150.5 && c<250.5)
        AQI=Linear(300,201,250.4,150.5,c)
    elseif (c>=250.5 && c<350.5)
        AQI=Linear(400,301,350.4,250.5,c)
    elseif (c>=350.5 && c<500.5)
        AQI=Linear(500,401,500.4,350.5,c)
    else
        AQI=501
    end
    return AQI
end

function AQIPM10(Concentration)
    Conc = Concentration
    ismissing(Conc) && return missing
    c = round(Int, Conc)
    if (c>=0 && c<55)
        AQI=Linear(50,0,54,0,c)
    elseif (c>=55 && c<155)
        AQI=Linear(100,51,154,55,c)
    elseif (c>=155 && c<255)
        AQI=Linear(150,101,254,155,c)
    elseif (c>=255 && c<355)
        AQI=Linear(200,151,354,255,c)
    elseif (c>=355 && c<425)
        AQI=Linear(300,201,424,355,c)
    elseif (c>=425 && c<505)
        AQI=Linear(400,301,504,425,c)
    elseif (c>=505 && c<605)
        AQI=Linear(500,401,604,505,c)
    else
        AQI=501
    end
    return AQI
end

function AQICO(Concentration)
    Conc = Concentration*0.873e-3
    ismissing(Conc) && return missing
    c = round(Conc; digits = 1)
    if (c>=0 && c<4.5)
        AQI=Linear(50,0,4.4,0,c)
    elseif (c>=4.5 && c<9.5)
        AQI=Linear(100,51,9.4,4.5,c)
    elseif (c>=9.5 && c<12.5)
        AQI=Linear(150,101,12.4,9.5,c)
    elseif (c>=12.5 && c<15.5)
        AQI=Linear(200,151,15.4,12.5,c)
    elseif (c>=15.5 && c<30.5)
        AQI=Linear(300,201,30.4,15.5,c)
    elseif (c>=30.5 && c<40.5)
        AQI=Linear(400,301,40.4,30.5,c)
    elseif (c>=40.5 && c<50.5)
        AQI=Linear(500,401,50.4,40.5,c)
    else
        AQI=501
    end
    return AQI
end

function AQISO21hr(Concentration)
    Conc=Concentration*0.382
    ismissing(Conc) && return missing
    try
    c=round(Int, Conc)
    if (c>=0 && c<36)
        AQI=Linear(50,0,35,0,c)
    elseif (c>=36 && c<76)
        AQI=Linear(100,51,75,36,c)
    elseif (c>=76 && c<186)
        AQI=Linear(150,101,185,76,c)
    elseif (c>=186 && c<=304)
        AQI=Linear(200,151,304,186,c)
    elseif (c>=304 && c<=604)
        AQI=201
    else
        AQI=501
    end
    return AQI
    catch e
        println("Error occured $e")
        println("Conc")
    end
end

function AQISO224hr(Concentration)
    Conc=Concentration*0.382
    ismissing(Conc) && return missing
    c=round(Int, Conc)
    if (c>=0 && c<=304)
        AQI=0
    elseif (c>=304 && c<605)
        AQI=Linear(300,201,604,305,c)
    elseif (c>=605 && c<805)
        AQI=Linear(400,301,804,605,c)
    elseif (c>=805 && c<=1004)
        AQI=Linear(500,401,1004,805,c)
    else
        AQI=501
    end
    return AQI
end

function AQIOzone8hr(Concentration)
    Conc=Concentration*0.509e-3
    ismissing(Conc) && return missing
    c=round(Conc; digits=3)
    if (c>=0 && c<.055)
        AQI=Linear(50,0,0.054,0,c)
    elseif (c>=.055 && c<.071)
        AQI=Linear(100,51,.070,.055,c)
    elseif (c>=.071 && c<.086)
        AQI=Linear(150,101,.085,.071,c)
    elseif (c>=.086 && c<.106)
    AQI=Linear(200,151,.105,.086,c)
    elseif (c>=.106 && c<.201)
        AQI=Linear(300,201,.200,.106,c)
    elseif (c>=.201 && c<.605)
        AQI=301
    else
        AQI=501
    end
    return AQI
end

function AQIOzone1hr(Concentration)
    Conc=Concentration*0.509e-3
    ismissing(Conc) && return missing
    c=round(Conc; digits=3)
    if (c>=0 && c<=.124)
        AQI=0
    elseif (c>=.125 && c<.165)
        AQI=Linear(150,101,.164,.125,c)
    elseif (c>=.165 && c<.205)
        AQI=Linear(200,151,.204,.165,c)
    elseif (c>=.205 && c<.405)
        AQI=Linear(300,201,.404,.205,c)
    elseif (c>=.405 && c<.505)
        AQI=Linear(400,301,.504,.405,c)
    elseif (c>=.505 && c<.605)
        AQI=Linear(500,401,.604,.505,c)
    else
        AQI=501
    end
    return AQI
end

function AQINO2(Concentration)
    Conc=Concentration*0.531e-3
    ismissing(Conc) && return missing
    c=round(Conc; digits=3)
    if (c>=0 && c<0.054)
        AQI=Linear(50,0,0.053,0,c)
    elseif (c>=0.054 && c<.101)
        AQI=Linear(100,51,.100,0.054,c)
    elseif (c>=.101 && c<.361)
        AQI=Linear(150,101,.360,.101,c)
    elseif (c>=.361 && c<.650)
        AQI=Linear(200,151,.649,.361,c)
    elseif (c>=.650 && c<1.250)
        AQI=Linear(300,201,1.249,.650,c)
    elseif (c>=1.250 && c<1.650)
        AQI=Linear(400,301,1.649,1.250,c)
    elseif (c>=1.650 && c<=2.049)
        AQI=Linear(500,401,2.049,1.650,c)
    else
        AQI=501
    end
    return AQI
end


"""

    c2a(conc, Pollutant)

Take input concentration in μg/m³ and the name of pollutant. Pollutant handled are "PM2.5", "PM10", "NO2", "CO", "SO2", "O3", "O38hr", "SO224hr"
Returns: AQI or missing
"""
function c2a(conc, Pollutant)
    func_dict = Dict("PM2.5"=>AQIPM25, "PM10"=>AQIPM10,
                    "NO2"=>AQINO2, "CO"=>AQICO, "SO2"=>AQISO21hr,
                    "O3"=>AQIOzone1hr, "O38hr"=>AQIOzone8hr,
                    "SO224hr"=>AQISO224hr)
    return func_dict[Pollutant](conc)
end

"""

    a2c(AQI, Pollutant)

Take input AQI and the name of pollutant. Pollutant handled are "PM2.5", "PM10", "NO2", "CO", "SO2", "O3", "O38hr", "SO224hr"
Returns: Conc in μg/m³ or missing
"""
function a2c(AQI, Pollutant)
    func_dict = Dict("PM2.5"=>ConcPM25, "PM10"=>ConcPM10,
                    "NO2"=>ConcNO2, "CO"=>ConcCO, "SO2"=>ConcSO21hr,
                    "O3"=>ConcOzone1hr, "O38hr"=>ConcOzone8hr,
                    "SO224hr"=>ConcSO224hr)
    return func_dict[Pollutant](AQI)
end

function ci2cf(Conc, Pollutant; AQI_diff = 10)
    ismissing(Conc) && return missing
    AQI_In = c2a(Conc, Pollutant)
    AQI_Out = AQI_In ≤ 50 ? min(50, AQI_In + AQI_diff) : AQI_In
#     AQI_Out = max(50, AQI_In)
    Conc_Out = a2c(AQI_Out, Pollutant)
    return max(Conc, Conc_Out)
end

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 16:07:06 2017

@author: jeremyinman
"""

import urllib
import csv
import json
import numpy as np
import math
from scipy.stats import poisson

def getdata(ticker):
    try:
        #store our dates
        datelist=[]
    
        #api key and url build
        key='TJLL4VNF9L2C2CNT'
        url="https://www.alphavantage.co/query?function=TIME_SERIES_DAILY_ADJUSTED&symbol="+ticker+"&outputsize=full&apikey="+key
        u=urllib.request.urlopen(url)
        content=u.read()
        
        #pull json data
        data=json.loads(content)
        timedata=data['Time Series (Daily)']
        #take all TS data
        
        #for the array size
        days=len(timedata)
        
        #using numpy will make calculations later faster
        closingprices=np.zeros(days,dtype=float)
        dividends=np.zeros(days,dtype=float)
        
        #map placement in array
        i=0
        
        #for each date in our TS, pull close, dividend and date
        for element in data['Time Series (Daily)']:
            x=timedata[element]
            datelist.append(element)
            close=float(x['5. adjusted close'])
            div=float(x['7. dividend amount'])
            closingprices[i]=close               
            dividends[i]=div
            i+=1
    
        return(datelist,closingprices,dividends)
    
    #handle failed requests... sometimes service fails to return data
    except urllib.error.HTTPError:
        
        datelist='lookup error'
        closingprices='lookup error'
        dividends='lookup error'
        
        return(datelist,closingprices,dividends)
    
    except KeyError:
    
        datelist='lookup error'
        closingprices='lookup error'
        dividends='lookup error'
        
        return(datelist,closingprices,dividends)
        
    
def calc_sigma(closingprices): 
    #need to calculate daily returns
    try:
        days=len(closingprices)
    
        returns=np.zeros(days,dtype=float)
    
        #counters for returns 
        j=0
        i=1
        for j in range(90):
            #also possible to use natural log
            returns[j]=(closingprices[i-1]/closingprices[i])-1
            i+=1
            j+=1
            
    except ValueError: #if div by zero or something weird, just give 0
        returns[j]=0
        i+=1
        j+=1
            
    except TypeError:
        returns=0
        
    except IndexError:
        returns=0
             
    try:
    #calculate recent volatilities 
        sigma15=np.std(returns[:15])/math.sqrt(15)
        sigma30=np.std(returns[:30])/math.sqrt(30)
        sigma60=np.std(returns[:60])/math.sqrt(60)
        sigma90=np.std(returns[:90])/math.sqrt(90)
        return(sigma15,sigma30,sigma60,sigma90)
    
    except TypeError:
        sigma15='type error'
        sigma30='type error'
        sigma60='type error'
        sigma90='type error'
        return(sigma15,sigma30,sigma60,sigma90)

def pois_jumps(closingprices):
    #need to calculate daily returns
    days=len(closingprices)
    returns=np.zeros(days,dtype=float)

    jumps_up=[]
    jumps_down=[]
    jump_up_criteria=0.1
    jump_down_criteria=-0.1
    
    #counters for returns 
    i=1
    j=0
    for j in range(days-1):
        try: #also possible to use natural log
            returns[j]=(closingprices[i-1]/closingprices[i])-1
            if returns[j]>=jump_up_criteria:
                jumps_up.append(returns[j])
                i+=1
                j+=1
            elif returns[j]<=jump_down_criteria:
                jumps_down.append(returns[j])
                i+=1
                j+=1
            else:
                i+=1
                j+=1
        except ValueError: #if div by zero or something weird, just give 0
            returns[j]=0
            j+=1
            
        except TypeError:
            returns[j]=0
            j+=1
    
    Jump_up=len(jumps_up)      
    Jump_down=len(jumps_down)
    jump_stats=[]
    jump_stats.append(np.mean(jumps_up))
    jump_stats.append(np.mean(jumps_down))
    jump_stats.append(np.std(jumps_up))
    jump_stats.append(np.std(jumps_down))
    return(Jump_up,Jump_down,jump_stats)

def calc_divYield(dividends,closingprices):
    try:
        #need number of dividends per year, curr price, and most recent dividend
        days=252
        current_price = closingprices[0]
    
        div_paid=[]
        i=0
        for i in range(days-1):
            if dividends[i]==0:
                i+=1
            else:
                div_paid.append(dividends[i])
                i+=1
    
        recent_div=div_paid[0]
    
        divs_per_year=len(div_paid)

        divYield = (divs_per_year*recent_div)/current_price
        return(divYield)
    
    except ValueError:
        divYield=0
        return(divYield)
    
    except IndexError:
        divYield=0
        return(divYield)
    
def monte_carlo_asset_valuation(closingprices,sigma,divYield,num_Jumps_up,num_Jumps_down,jump_stats,):
    
    try:
        n = 500 #time steps for each trading day
        Trials = 100 #simulations per run
        T=21/252 #time steps per day
        h= T/n
        
        r=.0272 #30 yr tbill aka 'risk free' rate
        
        days=len(closingprices)
        years=days/365
        months=years*12
        #Occurrences of positive asset price jumps, occurrences/years surveyed
        Lambda_up = num_Jumps_up/months
        Lambda_down = num_Jumps_down/months
        
        #Mean of jumps surveyed
        alpha_jump_up = jump_stats[0]
        alpha_jump_down = jump_stats[1]
        
        #Standard deviation of jumps surveyed
        sigma_jump_up = jump_stats[2]
        sigma_jump_down = jump_stats[3]
    
        #Poisson probability thresholds
        up_jump1_threshold = 1 - poisson.pmf(1,Lambda_up*h)
        up_jump2_threshold = 1 - poisson.pmf(2,Lambda_up*h)
        up_jump3_threshold = 1 - poisson.pmf(3,Lambda_up*h)
        up_jump4_threshold = 1 - poisson.pmf(4,Lambda_up*h)
        down_jump1_threshold = 1 - poisson.pmf(1,Lambda_down*h)
        down_jump2_threshold = 1 - poisson.pmf(2,Lambda_down*h)
        down_jump3_threshold = 1 - poisson.pmf(3,Lambda_down*h)
        down_jump4_threshold = 1 - poisson.pmf(4,Lambda_down*h)

        Expiration_Assets = np.zeros(Trials,dtype=np.float32)
        rands=np.zeros((Trials+1,n+1),dtype=np.float32)
    
        trial=1
        Jumps_Up=0
        Jumps_Down=0
    
        while trial < Trials + 1: 
         
            S=np.zeros(n+1,dtype=np.float32)
            N=np.zeros(n+1,dtype=np.float32)        
        
            i=0 #time step counter
    
            for i in range(0,n+1):
                N[i] = i*(T/n)     
        
            S[0] = closingprices[0]
    
            k=1 #time step counter
    
            for k in range(1,n+1):                
                
                m_up=0
                    
                draw=np.random.random()
                
                if draw > up_jump4_threshold:
                    m_up=4
                elif draw > up_jump3_threshold:
                    m_up=3
                elif draw > up_jump2_threshold:
                    m_up=2
                elif draw > up_jump1_threshold:
                    m_up=1
                else:
                    m_up=0
        
                Jumps_Up += m_up #keeping a running total of jumps up in all trials
                    
                m_down=0
                    
                draw=np.random.random()
            
                if draw > down_jump4_threshold:
                    m_down=4
                elif draw > down_jump3_threshold:
                    m_down=3
                elif draw > down_jump2_threshold:
                    m_down=2
                elif draw > down_jump1_threshold:
                    m_down=1
                else:
                    m_down=0
            
                Jumps_Down += m_down #keeping a running total of jumps down in all trials
                    
                #calcualte jump values 
                jumpAdjustment_up = Lambda_up*(math.exp(alpha_jump_up) - 1)
                jumpAdjustment_down = Lambda_down*(math.exp(alpha_jump_down) - 1)
            
                Jump_Drift_Up = math.exp(m_up*(alpha_jump_up - 0.5*sigma_jump_up**2))
                Jump_Drift_Down = math.exp(m_down*(alpha_jump_down - 0.5*sigma_jump_down**2))
                
                Jump_Noise_Rand_Sum = 0
                f=1
                
                while f < m_up+1:
                    Jump_Noise_Rand_Sum += np.random.normal(0,1)
                    f += 1
                
                Jump_Noise_Up = math.exp(sigma_jump_up*Jump_Noise_Rand_Sum)
                
                Jump_Noise_Rand_Sum = 0
                f=1
                
                while  f < m_down+1:
                    Jump_Noise_Rand_Sum += np.random.normal(0,1)
                    f += 1
                
                Jump_Noise_Down = math.exp(sigma_jump_down*Jump_Noise_Rand_Sum)                                
                   
                #anti-thetic variate
                #in order to prevent additional drift, basically
                #capture the first half of the random numbers
                #and use the negative of them for the second half
                if trial < Trials/2+1:
                    rands[trial][k]=np.random.normal(0,1)
                else:
                    rands[trial][k]=-1*rands[trial-int(Trials/2)][k]
                    
                #calculate our drift and noise terms
                Drift = math.exp((r - jumpAdjustment_up - jumpAdjustment_down - divYield - 0.5*sigma**2)*h)
                Noise = math.exp(sigma*math.sqrt(h)*rands[trial][k])
                
                #implement the asset price calculator for each time step
                S[k] = S[k-1]*Drift*Noise*Jump_Drift_Up*Jump_Noise_Up*Jump_Drift_Down*Jump_Noise_Down
                 
                     #capture end price of each simulation
                if k == n:
                    Expiration_Assets[trial-1] = S[k]
                    trial+=1
                    
        #calculate price
        Average_Expiration_Asset_Price = np.mean(Expiration_Assets)
        return(Average_Expiration_Asset_Price)
    
    #possible errors
    except ValueError:
        
        Average_Expiration_Asset_Price = "###"
        return(Average_Expiration_Asset_Price)
    
    except UnboundLocalError:
        
        Average_Expiration_Asset_Price = "###"
        return(Average_Expiration_Asset_Price)
    
    except TypeError:
        Average_Expiration_Asset_Price = "###"
        return(Average_Expiration_Asset_Price)        

def tickerlists(tickerlist):
    
    with open(tickerlist,'r') as f:
        Tickers=[x.strip() for x in f.readlines()]
        return Tickers
    f.close()
    
#file name
fname="valuations.dat" #output file name
#choose your exchange or other list from available options below
f_NASDAQ="nasdaq_symbols.txt"
f_NYSE="nyse_symbols.txt"
f_AMEX="amex_symbols.txt"
f_sp500="sp500.txt"

tickers=["NASDAQ:AAPL"]
Tickers=tickerlists(f_sp500) #enter variable name for list here
numb=len(Tickers)

#start program here
i=0 
with open(fname,'a') as f: #create our data file
       
    writer=csv.writer(f,dialect="excel")
    
    for i in range(len(Tickers)): #for each stock ticker
        try:  
                
            print(Tickers[i], " Loading...") #print outs to know its working
            data=getdata(Tickers[i]) #retrieve the data via api
        
            dates=data[0] #sort the data into our necessary parts 
            closing=data[1]
            divs=data[2]
        
            sigmas=calc_sigma(closing) #get our volatility calculations - we have 15,30,60, and 90 day vols
        
            poisson_jumps=pois_jumps(closing) #get our pois jump data
        
            divYield=calc_divYield(divs,closing) #calc divyield
        
            j=0 #counter for simulator
        #inputs for monte: closing prices, vol, divyield, pois jumps[0,1,2]
            Asset_Value=monte_carlo_asset_valuation(closing,sigmas[3],divYield,poisson_jumps[0],poisson_jumps[1],poisson_jumps[2])
        
        #return a value
            print(test, '--',Tickers[i], '--', Asset_Value, '--', closing[0], '--', 'Change: ', Asset_Value-closing[0])
        
            #dictate what we want to dump in the .dat file
            output=(Tickers[i],Asset_Value,dates,closing,divs,sigmas,poisson_jumps,divYield)
            writer.writerow(output) #write the data
            test+=1
        
        except TypeError:
            print(test, '--', Tickers[i], 'Type Error')
            output=(Tickers[i],Asset_Value,dates,closing,divs,sigmas,poisson_jumps,divYield, 'type error' )
            writer.writerow(output)
            test+=1
f.close()

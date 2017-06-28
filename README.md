# Mass-Asset-Valuation
Valuing companies on a large scale can be difficult. Luckily, you can use data from daily returns and dividends in order to run an options-based valuation. This program is designed to use an options-based valuation approach and implement it for any symbol in a list (i.e. an exchange or the S&amp;P 500). 

Using Alpha Vantage's free api, you can retrieve daily close prices for any company. So simply create a list of tickers to run through and you can request data. Alpha Vantage provides documentation on the data output and request formatting on their site. Your data will be in JSON. Using a JSON reader, you can extract the daily closes, dates, and dividends from them.

With this data, you can do a lot! I choose to run it through a Monte-Carlo based Option-valuation algorithm based off of Black-Scholes models. Using this to get an estimated asset price for each firm at a time period of your choosing.

Use the Value_Whole_Exchange.py file to run this. I have provided the nasdaq and SP500 tickers as well. You can also see the montecarlo model as well. I run these through a Python Console in Anacond's Spyder editor. The data output file can be set but I have it default to valuation.dat 

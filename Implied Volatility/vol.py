import pandas as pd
import numpy as np
import scipy.stats as scs


class BlackScholes:
    '''
    This is a class for Options contract for pricing European options on stocks without dividends.
    
    Attributes:
        spot                    - required   : spot price of underlying asset (int or float)
        strike                  - required   : strike of the option (int or float)
        rate                    - required   : risk free rate (float)
        days_to_maturity        - required   : days to expiration (int or float)
        volatility              - required   : volatility of underlyring asset (float)        
    '''
    
    def __init__(self, spot, strike, rate, days_to_maturity, volatility,callprice =None, putprice =None):
        self.spot = spot
        self.strike = strike
        self.rate = rate
        self.days_to_maturity = days_to_maturity
        self.volatility = volatility
        self.callprice = callprice                 # market price of call
        self.putprice = putprice                   # market price of put
        
        # other instance variables, i have used underscore just for clarity - not necessary
        
        self._a_  = self.volatility * np.sqrt(self.days_to_maturity)                        # sigma * sqrt(t)
        
        if self.strike ==0:
            raise ZeroDivisionError                                                         # ZeroDivisonError
        else:
            self._d1_ = (1/self._a_) * (np.log(self.spot / self.strike) + \
                        (self.rate + (self.volatility**2) / 2) * self.days_to_maturity)
        
        self._d2_ = self._d1_ - self._a_
        self._b_ = np.exp(-self.rate * self.days_to_maturity)                               # exp(-rt)
        
    
    
        # storing all the values in list using private methods
        [self.call_price, self.put_price] = self._price()
        [self.call_delta, self.put_delta] = self._delta()
        self.gamma = self._gamma()
        self.vega = self._vega()
        self.imp_vol = self._impvol()                                                       # implied volatility
        [self.call_theta, self.put_theta] = self._theta()
        [self.call_rho, self.put_rho] = self._rho()
    
    
    
    # underscore is used to make this function as private and not  public, can't call this using class object
        
    # option price
    def _price(self):
        '''Returns the option price: [Call price, Put price]'''

        # if volatility is zero or days to expiry is zero means option is at maturity
        if (self.volatility == 0) | (self.days_to_maturity == 0):
            call = max(self.spot - self.strike,0.0)
            put = max(self.strike - self.spot,0.0)
        else:
            call_price = self.spot * scs.norm.cdf(self._d1_) - self.strike * self._b_ * scs.norm.cdf(self._d2_)
            put_price = self.strike * self._b_ * scs.norm.cdf(-self._d2_) - self.spot * scs.norm.cdf(-self._d1_)

        return [call_price,put_price]

    # option delta
    def _delta(self):
        '''Returns the option delta: [Call delta, Put delta]'''

        # if volatility is zero or days to expiry is zero means option is at maturity
        if (self.volatility == 0) | (self.days_to_maturity == 0):
            call_delta = 1.0 if self.spot > self.strike else 0.0
            put_delta = -1.0 if self.strike > self.spot else 0.0
        else:
            call_delta = scs.norm.cdf(self._d1_)
            put_delta = -scs.norm.cdf(-self._d1_)

        return [call_delta,put_delta]

    # option gamma
    def _gamma(self):                                                                                    
        '''Returns the option gamma'''

        gamma = scs.norm.pdf(self._d1_)/(self.spot*self._a_)

        return gamma

    # option vega
    def _vega(self):
        '''Returns the option vega'''

        vega = self.spot*scs.norm.pdf(self._d1_)*np.sqrt(self.days_to_maturity)

        return vega

    # option theta
    def _theta(self):               
        '''Returns the option theta [Call theta, Put theta]'''

        call_theta = (-self.spot*scs.norm.pdf(self._d1_)*self.volatility/(2*np.sqrt(self.days_to_maturity))
                      -self.rate*self.strike*self._b_*scs.norm.cdf(self._d2_))
        put_theta = (-self.spot*scs.norm.pdf(self._d1_)*self.volatility/(2*np.sqrt(self.days_to_maturity))
                     + self.rate*self.strike*self._b_*scs.norm.cdf(-self._d2_))

        return [call_theta, put_theta]          

    # option rho
    def _rho(self):
        '''Returns the option rho [Call rho, Put rho]'''

        call_rho = self.strike*self.days_to_maturity*self._b_*scs.norm.cdf(self._d2_)
        put_rho = -self.strike*self.days_to_maturity*self._b_*scs.norm.cdf(-self._d2_)

        return [call_rho, put_rho]
    
    
    # option implied volatility
    def _impvol(self):
        '''Returns the option implied volatility '''
        
        if (self.callprice or self.putprice) is None:
            return self.volatility
        else:
            def bisection(high =500, low= 0):

                tolerance = 1e-7
                max_iter = 200

                if self.callprice:
                    price = self.callprice
                if self.putprice and not callprice:
                    price = self.putprice

                for i in range(max_iter):
                    mid = (high + low) / 2

                    # we check condition : implied volatility can't go negative

                    if mid < tolerance:
                        mid = tolerance

                    if self.callprice:
                        estimate = BlackScholes(self.spot,self.strike,self.rate,self.days_to_maturity,mid).call_price
                    if self.putprice:
                        estimate = BlackScholes(self.spot,self.strike,self.rate,self.days_to_maturity,mid).put_price

                    if round(estimate,6)== price:
                        break
                    elif estimate > price:
                        high = mid
                    elif estimate < price:
                        low = mid

                return mid
            
            return bisection()
                               

    
    
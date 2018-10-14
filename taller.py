# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 11:10:47 2018

@author: Juan
"""
import math as m
from scipy import optimize

class Derivada:
    def __init__(self, f, metodo="adelante", dx=0.001):
        self.f=f
        self.metodo=metodo
        self.dx=dx
        
    def adelante(self,x0):
        der=(self.f(x0+self.dx)-self.f(x0))/self.dx
        return der
    
    def central(self,x0):
        der=(self.f(x0+(self.dx/2))-self.f(x0-(self.dx/2)))/self.dx
        return der
    
    def extrapolada(self,x0):
        der=(4*self.central(x0)-self.central(x0))/3 
        return der
    
    def calc(self,x):
        if self.metodo=="adelante":
            calc=self.adelante(x)
        elif self.metodo=="central":
            calc=self.central(x)
        elif self.metodo=="extrapolada":
            calc=self.extrapolada(x)
        return calc
    
class Zeros:
    def __init__(self, f, metodo, error=1e-4, max_iter=100):
        self.f=f
        self.metodo=metodo
        self.error=float(error)
        self.max_iter=max_iter
        self.der=Derivada(self.f, "extrapolada", dx=error)
        
    def newton(self, itera, vi):
        for i in range(itera):
            vi=vi-self.f(vi)/self.der.calc(vi)
            if m.fabs(self.f(vi))<self.error:
                break
        return vi
    
    def interpolacion(self, itera, vi):
        a=vi[0]
        b=vi[1]
        for i in range(itera):
            if self.f(b)>0:
                c=a
                a=b
                b=c
            x=((b-a)/(self.f(b)-self.f(a)))*(-self.f(a))+a
            if self.f(x)>0:
                a=x
            elif self.f(x)<0:
                b=x
            else:
                break
            if m.fabs(self.f(a))<self.error or m.fabs(self.f(b))<self.error:
                break
            
        if m.fabs(self.f(a))>m.fabs(self.f(b)):
            return b
        else:
            return a
                
        
    def bisectriz(self, itera, vi):
        a=vi[0]
        b=vi[1]
        for i in range(itera):
            if self.f(b)>0:
                c=a
                a=b
                b=c
            x=(a+b)/2
            if self.f(x)>0:
                a=x
            elif self.f(x)<0:
                b=x
            else:
                break
            if m.fabs(self.f(a))<self.error or m.fabs(self.f(b))<self.error:
                break
        if m.fabs(self.f(a))>m.fabs(self.f(b)):
            return b
        else:
            return a
        
    def newton_sp(self, x0):
        res=optimize.newton(func=self.f, x0=x0, tol=self.error, maxiter=self.max_iter)
        return res
    
    def fsolve_sp(self, x0):
        res=optimize.fsolve(func=self.f, x0=x0, xtol=self.error, maxfev=self.max_iter)
        return res
    
    def brentq_sp(self, x0):
        res=optimize.brentq(f=self.f, a=x0[0], b=x0[1], xtol=self.error, maxiter=self.max_iter)
        return res
        
    def zero(self, vi):
        if type(vi)==int or type(vi)==float:
            vi=float(vi)
        elif type(vi)==tuple:
            vi=(float(vi[0]),float(vi[1]))
        if self.metodo=="newton":
            zero=self.newton(self.max_iter, vi)
        elif self.metodo=="interpolacion":
            zero=self.interpolacion(self.max_iter, vi)
        elif self.metodo=="bisectriz":
            zero=self.bisectriz(self.max_iter, vi)
        elif self.metodo=="newton-sp":
            zero=self.newton_sp(vi)
        elif self.metodo=="fsolve-sp":
            zero=self.fsolve_sp(vi)
        elif self.metodo=="brentq-sp":
            zero=self.brentq_sp(vi)
        return zero

if __name__=="__main__":
    x=Derivada(m.exp, "extrapolada")
    print(x.calc(1))
    y=Zeros(m.cos, "bisectriz", error=0.000001)
    print(y.zero((3,5)))
        
        
        
        
        
        
        

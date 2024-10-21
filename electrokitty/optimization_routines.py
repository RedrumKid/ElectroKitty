# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 10:28:55 2024

@author: ozbejv
"""

import numpy as np
import cma

class electrokitty_optimization_algorithms():
    """
    Class containing optimisation routines that are used when fitting parameters to data
    
    Implemented Nelder-Mead and CMA-ES from the cma package
    
    CMA-ES-inbuilt is available, but should be used for educational purposes only, as the cma
    implementaion offers more robust estimation
    """
    
    def nelder_mead(self, f, x_start,
                    step=0.1, no_improve_thershold=10**-9, 
                    no_improv_break=20, max_iter=10**3,
                    alpha=1, gamma=2, rho=0.5, sigma=0.5):
        # x_start=np.array(x_start)
        
        dim=len(x_start)
        prev_best=f(x_start)
        no_improv=0
        res=[[x_start,prev_best]]
        
        print("Started Optimization via Nelder-Mead")
        print()
        for i in range(dim):
            x=x_start.copy()
            if x[i] == 0.0:
                x[i]=10**-7
            x[i] += step*x[i]
            score=f(x)
            res.append([x,score])

        
        iters=0
        while 1:
            
            res.sort(key=lambda x: x[1])
            best=res[0][1]
            
            if max_iter and iters >= max_iter:
                print(f"fun evals: {iters*4}   current best: {res[0][1]}")
                print(f"current best: {res[0][0]}")
                print()
                print("Reached max num of iterations")
                return res[0], "Reached max num of iterations"
            
            iters+=1

            if best < prev_best-no_improve_thershold:
                no_improv=0
                prev_best=best
            else:
                no_improv+=1
            
            if no_improv >= no_improv_break:
                print(f"fun evals: {iters*4}   current best: {res[0][1]}")
                print(f"current best: {res[0][0]}")
                print()
                print("No improvevment after "+str(no_improv)+" iterations")
                return res[0], "No improvevment after "+str(no_improv)+" iterations"
            
            x0=[0]*dim
            
            for tup in res[:-1]:
                for i,c in enumerate(tup[0]):
                    x0[i] += c/(len(res)-1)

            xr=x0+alpha*(x0-res[-1][0])
            rscore=f(xr)
            
            if res[0][1] <= rscore < res[-2][1]:
                res[-1]=[xr,rscore]
                continue
            
            if rscore < res[0][1]:
                xe = x0 +gamma*(x0-res[-1][0])
                escore = f(xe)
                if escore < rscore:
                    res[-1]=[xe,escore]
                    continue
                else:
                    res[-1]=[xr, rscore]
                    continue
            
            xc = x0 + rho*(x0-res[-1][0])
            cscore = f(xc)
            if cscore < res[-1][1]:
                res[-1]=[xc, cscore]
                continue
            
            x1 = res[0][0]
            nres=[]
            for tup in res:
                redx = x1 + sigma*(tup[0]-x1)
                score=f(redx)
                nres.append([redx, score])
            
            res = nres
            print(f"fun evals: {iters*4}   current best: {res[0][1]}")
            print(f"current best: {res[0][0]}")
            print()


    def evaluate_func(self, check, lb, ub, xfeasable, x, fun, fmax):
        if check:
            if np.all(x>lb) and np.all(x<ub):
                return fun(x)
            else:
                return fmax+np.linalg.norm(x-xfeasable)
        else:
            return fun(x)
        
    def check_for_stop(self,check_arguments, mean, cov_mat, pc, sigma, numfmin, D, fbest, iteration):
        message="No improv after full iteration scope"
        stop=False
        success=False
        
        # NoEffectCoord
        if np.any(mean==mean+0.2*sigma*np.diag(cov_mat)):
            stop=True
            message="No improv in mean update"
            success=False
            return stop, message, success
        
        elif numfmin >= check_arguments[0]:
            stop=True
            message=f"No improvment on optimization in {numfmin} evaluations"
            success=False
            return stop, message, success
        
        elif np.max(sigma)*np.max(D) - check_arguments[1] > 10**4:
            stop=True
            message="sigma*max(diag(D)) increased by 10^4"
            success=False
            return stop, message, success
        
        elif iteration >= check_arguments[2] and fbest <= check_arguments[3]:
            stop =True
            message="f_best under tolerance"
            success=True
            return stop, message, success
        
        elif np.all(np.diag(cov_mat)<check_arguments[4]) and np.all(sigma*pc < check_arguments[4]):
            stop=True
            message="covariance close to delta distribution"
            success=True
            return stop, message, success
        else:
            return stop, message, success
        
    def single_CMA_ES_eval(self, fun, x, sigma,lb=None, ub=None, tolf=10**-12, tolc=10**-12, print_val=10):
        
        bounded_optimization=False
        xfeasable=x
        if lb !=None and ub != None:
            lb=np.array(lb)
            ub=np.array(ub)
            bounded_optimization=True
            xfeasable=(ub-lb)/2

        N=len(x)
        xmean=x
        fmax=fun(x)
        stopeval=10**3*N**2
        
        lambd=4+int(3*(np.log10(N)+1))
        mu=lambd/2
        weights=np.log(mu+0.5) - np.log(np.linspace(1,mu, int(mu)))
        weights=weights/(np.sum(weights))
        mueff=np.sum(weights)**2/np.sum(weights**2)
        mu=int(mu)
        
        cc=(4+mueff/N)/(N+4+2*mueff/N)
        cs=(mueff+2)/(N+mueff+5)
        c1=2/((N+1.3)**2+mueff)
        cmu=min([1-c1, 2*(mueff-2+1/mueff)/((N+2)**2+2*mueff/2)]) #+1/4
        damps=1+2*max([0, np.sqrt((mueff-1)/(N+1))-1])+cs
        
        pc=np.zeros(N)
        ps=np.zeros(N)
        B=np.eye(N)
        D=np.eye(N)
        C=np.eye(N)
        eigenval=0
        chiN=N**0.5*(1-1/(4*N)+1/(21*N**2))
        
        arx=np.zeros((N,lambd))
        arz=np.zeros((N,lambd))
        arfitness=np.ones(lambd)
        
        countereval=0
        best_solution=[x, fmax, countereval]
        check_arguments=[10+np.ceil(30*N/lambd), 1, 10+np.ceil(30*N/lambd), tolf, tolc*sigma]
        numfmin=0

        # stopeval=12
        
        while countereval < stopeval:
            xold=xmean
            C_old=C
            D_old=D
            B_old=B
        
            for k in range(lambd):
                arz[:,k]=np.random.multivariate_normal(np.zeros(N),np.eye(N))
                arx[:,k]=xmean + sigma*np.matmul(B,np.matmul(D,arz[:,k]))
                arfitness[k]=self.evaluate_func(bounded_optimization, lb, ub, xfeasable, 
                                           arx[:,k], fun, fmax)
                
                countereval+=1
            
            sort_ind=np.argsort(arfitness)
            
            if abs(best_solution[1] - arfitness[sort_ind[0]]) <= tolf:
                numfmin+=1
            
            if arfitness[sort_ind[0]] > fmax:
                fmax=arfitness[sort_ind[0]]
                
            xprop=0
            for i in range(mu):
                xprop+=+weights[i]*arx[:, sort_ind[i]]
            
            y=xmean-xold
            ps=(1-cs)*ps + np.sqrt(cs*(2-cs)*mueff)/sigma * np.matmul(B*D,y)
            if np.linalg.norm(ps)/np.sqrt(1-(1-cs)**(2*countereval/lambd))/chiN < 1.4+2/(N+1):
                hsig=1
            else:
                hsig=0
            
            pc=(1-cc)*pc + hsig*np.sqrt(cc*(2-cc)*mueff)/sigma*y
            
            C=(1-c1-cmu-c1*(1-hsig)*cc*(2-cc))*C + c1*(np.outer(pc,pc.T))  

            for i in range(mu):
                y=arx[:,sort_ind[i]]-xold
                C+=weights[i]*cmu/sigma**2*np.outer(y, y.T)


            sigma = sigma*min(1,np.exp((cs/damps)*(np.linalg.norm(ps)/chiN-1)))

            if countereval-eigenval > lambd/(c1+cmu)/N/10:
                eigenval=countereval
                C=np.triu(C)+np.triu(C, k=1).T
                D, B = np.linalg.eig(C)
                if np.any(D<0):
                    C=C_old
                    D=D_old
                    B=B_old
                else:
                    D=np.diag(np.sqrt(D))
            
            if best_solution[1] > arfitness[sort_ind[0]]:
                best_solution[0] = arx[:, sort_ind[0]]
                best_solution[1] = arfitness[sort_ind[0]]
                best_solution[2] = countereval
            
            xmean=xprop
                    
            check_arguments[1]=np.max(sigma)*np.max(D_old)
            stop, message, success = self.check_for_stop(check_arguments, 
                                                    xmean, C, pc, sigma, 
                                                    numfmin, D, best_solution[1], int(countereval/lambd))
            
            if stop:
                break
            
            if arfitness[sort_ind[0]] == arfitness[int(np.ceil(0.7*lambd))]:
                sigma = sigma*np.exp(0.2+cs/damps)
            
            print(f"fun evals: {int(countereval)}   current best: {arfitness[sort_ind[0]]}")
            print(f"current best: {best_solution[0]}")
            print(f"best fun: {best_solution[1]}")
            print()
            
        return best_solution, message, success

    def CMA_controller(self,fun, x, sigma, lb=None, ub=None, tolf=10**-15, tolx=10**-12, 
                       n_tries=5, restart=True):
        """
        Wrapper for CMA-ES-inbuilt
        """
        
        print("Started optimization via CMA-ES")
        print()
        trie=1
        print(f"Current try: {trie}")
        print("Optimizer working...")
        print()
        sol, message, success=self.single_CMA_ES_eval(fun, x, sigma, lb=lb, ub=ub, tolf=tolf, tolc=tolx)

        print()
        print(message)
        print()
        if restart:
            while not success:
                trie+=1
                print("starting anew")
                print(f"Current try: {trie}")
                print("Optimizer working...")
                print()
                sol, message, success=self.single_CMA_ES_eval(fun, np.array(sol[0]), 0.1**trie*np.ones(len(sol[0])), 
                                                              lb=lb, ub=ub,tolf=tolf, tolc=tolx)
                print(message)
                print()
                
                if trie >= n_tries:
                    message="Exceeded number of available tries"
                    print()
                    print("Error:")
                    print("Exceeded number of available tries")
                    print()
                    break
        
        return sol, message

    def CMA_ES(self, fun, x, sigma,
               lb=None, ub=None, tolf=10**-15, tolx=10**-12):
        """
        Wrapper for cma implemntation of CMA-ES
        """
        best_guess, es = cma.fmin2(fun, x, 0.1, 
                                   {"bounds": [lb, ub], 
                                    "tolfun": tolf, "tolx": tolx},
                                   )
        
        return [best_guess, es.result.fbest], "cma package succeded"
    
    
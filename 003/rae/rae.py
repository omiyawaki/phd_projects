import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
import misc.par as par

def rae(tau0, b, beta, Fs, Fa, plev, ps, n):
    
    tau = tau0*(plev/ps)**n
    
    Ts = ( ( ( Fs*(2+tau0) + Fa*(1+b*tau0/(b+1)) ) / (2+beta*tau0) ) / par.sigma )**(1/4)
    T = ( ( ((Fs+Fa-beta*par.sigma*Ts**4)*(1+tau)+Fa*b/tau0*(tau/tau0)**(b-1) - Fa*tau0/(b+1)*(tau/tau0)**(b+1)) / (2*(1-beta)) ) / par.sigma )**(1/4)

    OLR = ((Fs+Fa-beta*par.sigma*Ts**4)*(1) ) / (2*(1-beta))

    F1 = Fs+Fa*(1-(tau/tau0)**b)-beta*par.sigma*Ts**4
    
    if plev[1]-plev[0] > 0:
        Tlev0 = T[-1]
    else:
        Tlev0 = T[0]

    invstr = Tlev0 - Ts
    
    return T, Ts, invstr, OLR, F1

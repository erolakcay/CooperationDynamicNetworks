using StatsBase

### Function iterating the basic model, with fixed linking probabilities.

function networkiterate(netw, pn, pr, netsize)
    # one iteration of the basic model
    db = sample(1:netsize,2, replace=false); #sample a death and a birth, without replacement
    inherit = netw[:,db[2]].*(1.-(rand(netsize)-pn .> 0.)); # socially inherited connections
    randconn = (1-netw[:,db[2]]).*(1.-(rand(netsize)-pr .> 0.)); # random connections
    newconn = inherit+randconn #total connections
    newconn[db[2]]=1; #connect to the parent
    netw[:,db[1]]=newconn; #replace the dead individual with the newborn
    netw[db[1],:]=newconn; #same
    netw[db[1],db[1]]=0; #set self-link equal to zero
    return netw
end

function exppay(delta::Float64, pi::Array{Float64,1})
    (1+delta).^pi
end

function linpay(delta::Float64, pi::Array{Float64,1})
    clamp.(1 + delta*pi, 0, Inf)
end
### Function iterating the network over one time-step while playing the coauthor game, with fixed linking probabilities

function networkitCoauthor(netw::Array{Int64,2}, types::Array{Int64,1}, pn::Float64, pr::Float64, netsize::Int64, b::Float64, c::Float64, mu::Float64,d::Float64=0.0,delta::Float64=1.0,payfun::Function=linpay)
    death = rand(1:netsize); #sample a death at random
    degrees = vec(sum(netw,1)) #calculate degrees
    meandeg = mean(degrees)
    degrees[degrees.==0]=1 #set zero degree individuals to 1, to avoid division by zero (this has no effect to anyone's payoff, b/c of multiplication w/ adjacency matrix)
    outben = (types)./degrees #benefit going out from each individual
    inben = netw*outben #benefit flowing into each individual
    synben= inben.*outben #the product of in and out flows
    payoffs = payfun(delta, (b*inben - types*c + d*synben)) #calculate payoffs
    meanpayoff = mean(payoffs)
    payoffs[death]=0 #the dying individual gets 0 payoff (i.e., can't give birth)
    payoffw=weights(payoffs);
    birth = sample(1:netsize,payoffw) #sample parent with probability proportional to payoff
    ## The rest proceeds as before, except we need to assign type to the newborn
    newconn = netw[:,birth].*(1.-(rand(netsize)-pn .> 0.)) + (1-netw[:,birth]).*(1.-(rand(netsize)-pr .> 0.)) #total connections
    newconn[birth]=1; #connect to the parent
    netw[:,death]=newconn; #replace the dead individual with the newborn
    netw[death,:]=newconn; #same
    netw[death,death]=0; #set self-link equal to zero
    types[death] = mod(types[birth]+(rand()<mu),2) # copy the parent's type, mutate with probability mu
    return (netw,types, meandeg, meanpayoff)
end

### Function for one iteration of the network with evolving linking probabilities and playing the coauthor game, allows costly links

function networkitCoauthorEvolLink(netw::Array{Int64,2}, types::Array{Int64,1}, pns::Array{Float64,1}, prs::Array{Float64,1}, netsize::Int64, b::Float64, c::Float64, mu::Float64, mulink::Float64, sigmapn::Float64, sigmapr::Float64,clink::Float64=0.0, d::Float64=0.0, delta::Float64=1.0,payfun::Function=linpay)
    death = rand(1:netsize); #sample a death at random
    degrees = vec(sum(netw,1))
    meandeg = mean(degrees)
    degrees[degrees.==0]=1 #set zero degree individuals to 1, to avoid division by zero (this has no effect to anyone's payoff, b/c of multiplication w/ adjacency matrix)
    outben = (types)./degrees #proportional (to b) benefit going out from each individual
    inben = netw*outben
    synben= inben.*outben
    payoffs = payfun(delta, b*inben- types* c + synben*d - vec(sum(netw,1))*clink) #calculate payoffs
    meanpayoff = mean(payoffs)
    payoffs[death]=0 #the dying individual gets 0 payoff (i.e., can't give birth)
    payoffw=weights(payoffs)
    birth = sample(1:netsize,payoffw) #sample parent with probability proportional to payoff
    ## The rest proceeds as before, except we need to assign type to the newborn
    pns[death]=clamp.(pns[birth]+(rand()<mulink)*randn()*sigmapn,0,1) #copy parent's linking probabilities, with possible mutation
    prs[death]=clamp.(prs[birth]+(rand()<mulink)*randn()*sigmapr,0,1)
    newconn = netw[:,birth].*(1.-(rand(netsize)-pns[death] .> 0.)) + (1-netw[:,birth]).*(1.-(rand(netsize)-prs[death] .> 0.)) #total connections
    newconn[birth]=1; #connect to the parent
    netw[:,death]=newconn; #replace the dead individual with the newborn
    netw[death,:]=newconn; #same
    netw[death,death]=0; #set self-link equal to zero
    types[death] = mod(types[birth]+(rand()<mu),2) # copy the parent's type, mutate with probability mu
    return (netw,types,pns,prs, meandeg, meanpayoff)
end

## This function runs the Prisoner's dilemma game with constant per link benefits (as opposed to constant total benefits), but also linearly increasing costs per unit benefit given out.

function networkitPD(netw, types, pn, pr, netsize, b, c, mu, d=0.0,delta::Float64=1.0,payfun::Function=linpay)
    death = rand(1:netsize); #sample a death at random
    degrees = vec(sum(netw,1))
    meandeg = mean(degrees)
    inben = netw*types
    synben= inben.*types
    payoffs = payfun(delta, b*inben-(types*c).*degrees+d*synben) #calculate payoffs
    meanpayoff = mean(payoffs)
    payoffs[death]=0 #the dying individual gets 0 payoff (i.e., can't give birth)
    payoffw=weights(payoffs);
    birth = sample(1:netsize,payoffw) #sample parent with probability proportional to payoff
    ## The rest proceeds as before, except we need to assign type to the newborn
    newconn = netw[:,birth].*(1.-(rand(netsize)-pn .> 0.)) + (1-netw[:,birth]).*(1.-(rand(netsize)-pr .> 0.)) #total connections
    newconn[birth]=1; #connect to the parent
    netw[:,death]=newconn; #replace the dead individual with the newborn
    netw[death,:]=newconn; #same
    netw[death,death]=0; #set self-link equal to zero
    types[death] = mod(types[birth]+(rand()<mu),2) # copy the parent's type, mutate with probability mu
    return (netw,types, meandeg, meanpayoff)
end

### This is again the PD game, but with evolving linking probabilities

function networkitPDEvolLink(netw::Array{Int64,2}, types::Array{Int64,1}, pns::Array{Float64,1}, prs::Array{Float64,1}, netsize::Int64, b::Float64, c::Float64, mu::Float64, mulink::Float64, sigmapn::Float64, sigmapr::Float64,clink::Float64=0,d::Float64=0.0,delta::Float64=1.0,payfun::Function=linpay)
    death = rand(1:netsize); #sample a death at random
    degrees = vec(sum(netw,1))
    meandeg = mean(degrees)
    inben = netw*types
    synben= inben.*types
    payoffs = payfun(delta, b*inben-(types*c).*degrees+d*synben) #calculate payoffs
    meanpayoff = mean(payoffs)
    payoffs[death]=0 #the dying individual gets 0 payoff (i.e., can't give birth)
    payoffw=weights(payoffs)
    birth = sample(1:netsize,payoffw) #sample parent with probability proportional to payoff
    ## The rest proceeds as before, except we need to assign type to the newborn
    pns[death]=clamp.(pns[birth]+(rand()<mulink)*randn()*sigmapn,0,1) #copy parent's linking probabilities, with possible mutation
    prs[death]=clamp.(prs[birth]+(rand()<mulink)*randn()*sigmapr,0,1)
    newconn = netw[:,birth].*(1.-(rand(netsize)-pns[death] .> 0.)) + (1-netw[:,birth]).*(1.-(rand(netsize)-prs[death] .> 0.)) #total connections
    newconn[birth]=1; #connect to the parent
    netw[:,death]=newconn; #replace the dead individual with the newborn
    netw[death,:]=newconn; #same
    netw[death,death]=0; #set self-link equal to zero
    types[death] = mod(types[birth]+(rand()<mu),2) # copy the parent's type, mutate with probability mu
    return (netw,types,pns,prs, meandeg, meanpayoff)
end

### Function for running one run of the simulation, starting with a random network, with "burn-in" with fixed linking probabilities and then switching the coauthor game with or without evolving linking probabilities (as given by the funname)
### 07-03-17: added optional argument for changing how frequently the output is returned.
### 08-21-17: added optional argument for specifying whether the PD or coauthor game should be used in the evolving link or non-evolving link scenarios (default is coauthor)


function runSim(;pn::Float64=0.5, pr::Float64=0.1, netsize::Int64=100, generations::Int64=100, b::Float64=1.0, c::Float64=0.5, d::Float64=0.0, mu::Float64=0.01, evollink::Bool=false, mulink::Float64=0.0, sigmapn::Float64=0.05, sigmapr::Float64=0.01,clink::Float64=0.0,retint::Int64=0, funnoevollink::Function=networkitCoauthor, funevollink::Function=networkitCoauthorEvolLink,delta::Float64=1.0,payfun::Function=linpay)

## Burn-in period
    ## Initialize the network randomly
    netw=rand([0,1],(netsize,netsize));
    netw = netw .* transpose(netw);
    netw = (1-diagm([1 for i=1:netsize])).*netw;
    for i in 1:(netsize*20) ##First iterate the network for 20 generations (i.e., 20*netsize time steps) neutrally to get a steady state network
      netw = networkiterate(netw,pn,pr,netsize)
    end
    types = rand([0,1],netsize) #initialize types at random

if retint == 0
  retint = netsize
end


## Switch on the coauthor game
    ## If argument evollink is true, evaluate with evolving linking probabilities, if false, with fixed linking probabilities
    if evollink == true
      pnhist=zeros(generations); ##initialize arrays to record history
      prhist=zeros(generations);
      typehist = zeros(generations);
      degreehist = zeros(generations);
      payoffhist = zeros(generations);
      prs = clamp.(randn(netsize)*sigmapr+pr,0,1);
      pns = clamp.(randn(netsize)*sigmapn+pn,0,1);
      for i in 1:(generations*retint)
        netw, types, pns, prs, mdeg, mpay = funevollink(netw,types,pns,prs,netsize,b,c,mu,mulink,sigmapn,sigmapr,clink,d,delta,payfun)
        if mod(i,retint)==0
          pnhist[div(i,retint)]=mean(pns)
          prhist[div(i,retint)]=mean(prs)
          typehist[div(i,retint)]=mean(types)
          degreehist[div(i, retint)]=mdeg
          payoffhist[div(i, retint)]=mpay
        end
      end
      return (typehist, pnhist, prhist, degreehist, payoffhist)
    else
      typehist = zeros(generations); ##initialize arrays to record history
      degreehist = zeros(generations);
      payoffhist = zeros(generations);
      for i in 1:(generations*retint)
        netw, types, mdeg, mpay = funnoevollink(netw,types,pn,pr,netsize,b,c,mu,d,delta,payfun)
        if mod(i,retint)==0
          typehist[div(i,retint)]=mean(types)
          degreehist[div(i, retint)]=mdeg
          payoffhist[div(i, retint)]=mpay
        end
      end
      return (typehist, pn, pr, degreehist, payoffhist)
    end
end

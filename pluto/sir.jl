### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ f35b164e-fffd-11ed-07ed-736855b6df20
md"
# SIR model with Pluto.jl
## Simon D.W. Frost

Version 0.1.0 2023-05-31

https://github.com/epirecipes

This is a Pluto.jl port of the Rmarkdown of the SIR model, written by Ottar N. Bjørnstad.

The basic equations for the flow of hosts between **S**usceptible, **I**nfectious and **R**ecovered compartments are:
 
$\begin{aligned}
    \frac{dS}{dt} =& \underbrace{\mu N}_{\mbox{birth}} - \underbrace{\beta I \frac{S}{N}}_{\mbox{infection}} - \underbrace{\mu S}_{\mbox{death}} \label{eq:sirs}\\
     \frac{dI}{dt} =& \underbrace{\beta I \frac{S}{N}}_{\mbox{infection}} - \underbrace{\gamma I}_{\mbox{recovery}} - \underbrace{\mu I}_{\mbox{death}}  \label{eq:siri}\\
     \frac{dR}{dt} =& \underbrace{\gamma I}_{\mbox{recovery}} - \underbrace{\mu R}_{\mbox{death}} \label{eq:sirr}
\end{aligned}$

The assumptions of this version of the SIR model are:

- The infection circulates 
     in a population of size $N$, with a per capita background
     death rate, $\mu$, which is balanced by a
     birth rate $\mu N$. From the sum of the equations $dN/dt=0$ and $N=S+I+R$ is constant.

- The infection causes morbidity (not mortality).

- Newborns are recruited directly into the susceptible
     class at birth.

- Transmission of infection from infectious to
     susceptible individuals is controlled by a bilinear contact
     term $\beta I \frac{S}{N}$, from the assumption that the $I$ infectious individuals are independently and randomly
     mixing with all other individuals, so the fraction $S/N$ is with susceptible individuals; $\beta$ is the transmission rate.

- Infected individuals move directly into the 
     the infectious class and remains there 
     for an average infectious period of $1/(\gamma+\mu)$ time units.

- Recovered individuals are immune
     from re-infection for life.

For the basic SIR model $R_0 = \frac{\beta}{\gamma + \mu}$. The isoclines (sometimes called the null-clines)  are given by the solution to the 
equations $dS/dt=0$ and $dI/dt=0$ and partitions the phase plane into regions 
were $S$ and $I$ are increasing and decreasing. 
For $N=1$, the $I$-isocline is $S = (\gamma +\mu)/\beta = 1/R_0$
and the S-isocline is $I= \mu (1/S-1)/\beta$.

The resonant period is $\frac{2 \pi}{\mu (R_0 -1) (\mu + \gamma)}$. 
"

# ╔═╡ 8761c521-77b1-4009-8238-6d6ba3f82e0b
function sirmod(u, p, t)
    S,I,R = u
    β = p.β
    μ = p.μ
    γ = p.γ
    N = p.N
    dS = μ*(N-S) - β*S*I/N
    dI = β*S*I/N - (μ+γ)*I
    dR = γ*I - μ*R
    [dS, dI, dR]
end;

# ╔═╡ Cell order:
# ╟─f35b164e-fffd-11ed-07ed-736855b6df20
# ╟─8761c521-77b1-4009-8238-6d6ba3f82e0b

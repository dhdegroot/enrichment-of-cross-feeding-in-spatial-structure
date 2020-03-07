#### Population statistics

We describe a microbial community of three species: cooperator $A$, cooperator $B$, and cheater $C$. The fractions of these species in the initial population is denoted by $f_a, f_b, f_c,$ respectively. The cells in the population are distributed over droplets (micro-environments), where the number of cells in one droplet is poisson-distributed with mean $\lambda$:

$P(N_t = n_t) = \frac{\lambda^{n_t} e^{-\lambda}}{n_t!}$

The species are distributed over droplets (micro-environments) according to a multinomial distribution. So, given that a certain droplet has $n_t$ cells, we have

$P(N_a = n_a,N_b=n_b, N_c=n_t-n_a-n_b~|~N_t =n_t ) = \binom{n_t}{n_a,n_b,n_t-n_a-n_b} f_a^{n_a} f_b^{n_b} f_c^{n_t-n_a-n_b}$.

Therefore, the probability of finding a certain droplet is given by

$P(N_a = n_a,N_b=n_b, N_c=n_t-n_a-n_b) = P(N_t = n_t) P(N_a = n_a,N_b=n_b, N_c=n_t-n_a-n_b~|~N_t =n_t ) =$

​		$= \frac{\lambda^{n_t} e^{-\lambda}}{n_t!} \binom{n_t}{n_a,n_b,n_t-n_a-n_b} f_a^{n_a} f_b^{n_b} f_c^{n_t-n_a-n_b}$

Consequently, the average number of cells $n_a$ in a droplet is given by summing 

$\langle n_a\rangle = \sum_{n_t=1}^\infty \sum_{n_a = 0}^{n_t} \sum_{n_b = 0}^{n_t-n_a} n_a \frac{\lambda^{n_t} e^{-\lambda}}{n_t!} \binom{n_t}{n_a,n_b,n_t-n_a-n_b} f_a^{n_a} f_b^{n_b} f_c^{n_t-n_a-n_b}$

​		$= \sum_{n_t=1}^\infty  \frac{\lambda^{n_t} e^{-\lambda}}{n_t!} \left(\sum_{n_a = 0}^{n_t} \sum_{n_b = 0}^{n_t-n_a} n_a \binom{n_t}{n_a,n_b,n_t-n_a-n_b} f_a^{n_a} f_b^{n_b} f_c^{n_t-n_a-n_b}\right)$

​		$= \sum_{n_t=1}^\infty  \frac{\lambda^{n_t} e^{-\lambda}}{n_t!} n_t f_a$

​		$= f_a \sum_{n_t=1}^\infty n_t \frac{\lambda^{n_t} e^{-\lambda}}{n_t!}$

​		$= f_a \lambda$,

where we have used the average of the multinomial distribution from the 2nd to the 3rd line, and the average of the Poisson distribution from the 4th to the 5th line. Therefore, if we have $N_C$ droplets, we have in total $N_C f_a \lambda$ cells of species $A$. Now if we pick one cell of species $A$, what is the probability that it is in a droplet of a certain composition? Then, we have to know how many cells of species $A$ occur in such droplets, which is just the number of droplets of that type times the number of cells of species $A$ in that type of droplet. This is given by  

$N_C n_a \frac{\lambda^{n_t} e^{-\lambda}}{n_t!} \binom{n_t}{n_a,n_b,n_t-n_a-n_b} f_a^{n_a} f_b^{n_b} f_c^{n_t-n_a-n_b}$,

therefore the probability of a cell of species $A$ being in a specific droplet is (we denote this by $q^A_{n_a,n_b,n_c}$), is given by

$q^A_{n_a,n_b,n_c} = P(\text{cell of species A, being in droplet with }N_a = n_a,N_b = n_b,N_c = n_t -n_a-n_b)$,

​		$= \frac{1}{f_a \lambda} n_a \frac{\lambda^{n_t} e^{-\lambda}}{n_t!} \binom{n_t}{n_a,n_b,n_t-n_a-n_b} f_a^{n_a} f_b^{n_b} f_c^{n_t-n_a-n_b}$,

and similarly for species $B$ and $C$.

#### A formula for the average growth rate

In the droplets, the cells will start growing. Their eventual number of offspring depends on the initial frequencies of the different species. We denote the number of offspring $A$ made per initial cell of $A$ by $G^A$, which is thus dependent on the initial species frequencies $G^A= G^A(n_a,n_b,n_c)$. Furthermore, this is also dependent of the final total number of cells in the droplet, which we call the carrying capacity: $CC$. If, we have the formulae for $G^A, G^B, G^C$, we could calculate the average growth for each species by

$\langle G^A \rangle = \sum_{n_t=0}^\infty \sum_{n_a = 0}^{n_t} \sum_{n_b = 0}^{n_t-n_a} q^A_{n_a,n_b,n_c} G^A(n_a,n_b,n_t-n_a-n_b)$,

​		$=  \frac{1}{f_a \lambda} \sum_{n_t=0}^\infty  \frac{\lambda^{n_t} e^{-\lambda}}{n_t!} \sum_{n_a = 0}^{n_t} \sum_{n_b = 0}^{n_t-n_a} n_a \binom{n_t}{n_a,n_b,n_t-n_a-n_b} f_a^{n_a} f_b^{n_b} f_c^{n_t-n_a-n_b} G^A(n_a,n_b,n_t-n_a-n_b)$,

and similarly for species $B$ and $C$. 

#### Several alternatives for calculating the growth factors

To determine $G^A$ we make a few simplifying assumptions:

1. We assume that a certain part of the carrying capacity is facilitated by pairs of cooperators $A,B$. How large this carrying capacity is, is determined by the initial frequencies: $CC_{\text{comm}} = CC_0 \frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c}$. This means that, the more cheaters there are, the less efficient the available substrate is used, and less offspring can be made.
2.  The different species are differentially dependent of this facilitated carrying capacity. With a species-specific $CC^A_{\text{ind}}$, we can model that cells might even grow (although less) without pairs of cooperators around. Our experimental setting is best modelled by setting $CC^A_{\text{ind}} = CC^B_{\text{ind}}=0$, and $CC^{C}_{\text{ind}}>0$. In other words, the cooperators are compulsory cooperators and the cheater is an optional cheater.
3. The two cooperators are stoichiometrically coupled: they always end up in a 1-1 ratio.

Based on these assumptions, and with the aim of keeping the model as simple as possible, we ended up with.:

$G^A(n_a,n_b,n_c) = \frac{1}{n_a} \left(CC_0 \frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c} \right)\frac{1}{2} \left(1-\frac{n_c}{n_a+n_b+n_c}\right)$,

$G^B(n_a,n_b,n_c) = \frac{1}{n_b}\left(CC_0 \frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c} \right)\frac{1}{2} \left(1 -  \frac{n_c}{n_a+n_b+n_c}\right)$,

$G^C(n_a,n_b,n_c) = \frac{1}{n_c}\left(CC_0 \frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c} + CC^C_{\text{ind}}\right)\frac{n_c}{n_a+n_b+n_c}$.

**Another alternative**

$G^A(n_a,n_b,n_c) = \frac{1}{n_a} \left(CC_0 \frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c} \right) \frac{\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c}$,

$G^B(n_a,n_b,n_c) = \frac{1}{n_b}\left(CC_0 \frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c} \right)\frac{\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c}$,

$G^C(n_a,n_b,n_c) = \frac{1}{n_c}\left(CC_0 \frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c} + CC^C_{\text{ind}}\right)\left(1-\frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c}\right)$.

**Yet another alternative**

We can also take another approach. 

*Without cooperators*: In the absence of cooperator pairs, we assume that all species can grow to a cell number of $CC^{A/B/C}_{\text{ind}}$. 

$G^{A/B/C} (n_a,n_b,n_c) = \frac{1}{n_{a/b/c}} CC^{A/B/C}_\text{ind}$

*With lots of cooperators*: In the abundant presence of cooperator pairs, the cells together reach a certain carrying capacity $CC^0$. The eventual fraction of the cheater in this case is determined by linearly approximating the time evolution of a fraction in a competitive growth experiment $\dot{\phi}_i = k\phi_i(1-\phi_i)$, where we assume that the cheater grows faster than the cooperators. We approximate this by

 $\left(\phi_c\right)_\text{end} = \left(\phi_c\right)_\text{start} + k\left(\phi_c\right)_\text{start}\cdot\left(1-\left(\phi_c\right)_\text{start}\right)$,

where $\phi_c = \left(\frac{n_c}{n_a + n_b + n_c}\right)$. Then, we will determine $k$ by assuming that the maximal increase in this fraction should be about $10\%$, so that we need $k = 0.1/0.25 = 2/5$. So, we get

$G^{C}=\frac{1}{n_{c}}CC_0 \left((1+k)\left(\frac{n_c}{n_a + n_b + n_c}\right) - k\left(\frac{n_c}{n_a + n_b + n_c}\right)^2\right)$

$G^{A/B} = \frac{1}{n_{a/b}} CC_0 \frac{1}{2} \left(1-(1+k)\left(\frac{n_c}{n_a + n_b + n_c}\right) - k\left(\frac{n_c}{n_a + n_b + n_c}\right)^2\right)$

*With any number of cooperator pairs*: If we have some cooperator pairs, we will linearly interpolate between the above two extreme cases. The cooperator pair abundance is determined by $\frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c}$, since this captures the number of cells part of a cooperator pair compared to the number of cheaters. This abundance will determine our linear interpolation. We will thus get

$G^{C}=\left(1-\frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c}\right)\frac{1}{n_c} CC^{C}_\text{ind}+\frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c}\frac{1}{n_{c}}CC_0 \left((1+k)\left(\frac{n_c}{n_a + n_b + n_c}\right) - k\left(\frac{n_c}{n_a + n_b + n_c}\right)^2\right)$

$G^{A/B}=\left(1-\frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c}\right)\frac{1}{n_{a/b}} CC^{A/B}_\text{ind}+\frac{2\min\{n_a,n_b\}}{2 \min\{n_a,n_b\} + n_c}\frac{1}{n_{a/b}}CC_0\frac{1}{2} \left(1-(1+k)\left(\frac{n_c}{n_a + n_b + n_c}\right) - k\left(\frac{n_c}{n_a + n_b + n_c}\right)^2\right)$

which, admittedly, does not look pretty. However, if we interpret the terms it is not so bad

$G^{A/B/C} = \left(1-\text{coop abundance}\right)\cdot (\text{coop indep. growth}) + \left(\text{coop abundance}\right)\cdot(\text{coop dep. growth})$

#### Bringing it together

Plugging these equations for the growth factors in above, we can find the growth of the different species as a function of the parameters. We are mostly interested in the behaviour of $\langle G^A\rangle,\langle G^B\rangle,\langle G^C\rangle$ as a function of the average number of cells per droplet $\lambda$, and how this behaviour changes when we change the parameters $f_a,f_b,f_c$. 

#### Interpretation

- Cooperators have a fitness advantage at intermediate cells per compartment
- Check if this is robust to parameter changes. Is there a fixed lambda that seems to work?
- For this lambda, check if this still works after a few iterations with different starting fractions.


Bk - (Bk * sk * sk.T * Bk) / (sk.T * Bk * sk) + (yk * yk.T) / (yk.T * sk)

$$
\Delta {g^{(k)}} = {g^{(k + 1)}} - {g^{(k)}}\\
	{H_{k + 1}} = {H_k} + \frac{{\Delta {x^{\left( k \right)}}\Delta {x^{\left( k \right)T}}}}{{\Delta {x^{\left( k \right)T}}\Delta {g^{\left( k \right)}}}} - \frac{{\left[ {{H_k}\Delta {g^{\left( k \right)}}} \right]{{\left[ {{H_k}\Delta {g^{\left( k \right)}}} \right]}^T}}}{{\Delta {g^{\left( k \right)T}}{H_k}\Delta {g^{\left( k \right)}}}}
$$


$$
\Delta {g^{(k)}} = {g^{(k + 1)}} - {g^{(k)}}\\
	{H_{k + 1}} = {H_k} - \frac{{H_k}{\Delta {x^{\left( k \right)}}\Delta {x^{\left( k \right)T}}}{H_k}}{{\Delta {x^{\left( k \right)T}}{H_k}\Delta {x^{\left( k \right)}}}} +\frac{{\left[ {\Delta {g^{\left( k \right)}}} \right]{{\left[ {\Delta {g^{\left( k \right)}}} \right]}^T}}}{{\Delta {g^{\left( k \right)T}}\Delta {x^{\left( k \right)}}}}
$$

Hk = Hk + (sk-Hk*yk) * (sk-Hk*yk).T / ((sk-Hk*yk).T * yk)

$$
\Delta {g^{(k)}} = {g^{(k + 1)}} - {g^{(k)}}\\
	{H_{k + 1}} = {H_k} + \frac{(\Delta {x^{\left( k \right)}}-{H_k}\Delta {g^{\left( k \right)}})(\Delta {x^{\left( k \right)}}-{H_k}\Delta {g^{\left( k \right)}})^T} {(\Delta {x^{\left( k \right)}}-{H_k}\Delta {g^{\left( k \right)}})^T\Delta {g^{\left( k \right)}}}
$$

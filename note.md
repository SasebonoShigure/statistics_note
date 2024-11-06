# 数理统计笔记

## 第二章

### 基本分布

- 两点（Bernoulli）分布， $X_i\sim B(n,p)$
  - $P(X_i=x_i)=p^{x_i}(1-p)^{1-x_i}$
- 指数分布， $X_i\sim Exp(\lambda)$
  - $p(x)=\lambda e^{-\lambda x}$
  - $E(x)=\frac{1}{\lambda},Var(x)= \frac{1}{\lambda ^2}$
  - MLE: $\hat\lambda = \frac{1}{\bar x}$
- 正态分布， $X\sim N(\mu,\sigma^2)$
  - $p(x) = \frac{1}{\sqrt{2\pi \sigma^2}}e^{-\frac{(x-\mu)^2}{2\sigma^2}}$
  - $E(x)=\mu, Var(x)=\sigma^2$
  - MLE: $\hat\mu=\bar x,\hat\sigma^2=\frac{1}{n}\Sigma(x_i-\bar x)^2$
  - 样本方差 $S^2=\frac{1}{n-1}\Sigma(x_i-\bar x)^2$ 是 $\sigma^2$ 的最小方差无偏估计
  - 性质：
    - $X\sim N(\mu,\sigma^2),nX\sim N(n\mu,n^2\sigma^2)$
    - $X_1,\cdots,X_n\sim N(\mu,\sigma^2), i.i.d.,\Sigma X_i\sim N(n\mu,n\sigma^2),\bar X \sim N(\mu,\frac{\sigma^2}{n})$
- 均匀分布， $X\sim U[a,b]$
  - $E(x)=\frac{1}{b-a} , Var(x)=\frac{(b-a)^2}{12}$
  - MLE: $\hat a=\min\limits_{1\leq i\leq n} x_i,\hat b=\max\limits_{1\leq i\leq n} x_i$
  - 对于 $X_i \sim U[0,\theta]$ ，矩估计： $\tilde\theta=2\bar x$
- 卡方分布
  - $X_1,\cdots,X_n\sim N(0,1), i.i.d.,\Sigma_{i=1}^nX_i^2\sim\chi^2(n)$
  - $\chi^2(n)=\Gamma(\frac n 2,\frac1 2)$
  - 性质：
    - $\xi\sim\chi^2(m), \eta\sim\chi^2(n)$ ，二者独立， $\xi+\eta\sim\chi^2(m+n)$
    - $X\sim\chi^2(n),EX=n$
- Gamma分布
  - 与指数分布关系： $X\sim Exp(\lambda)=\Gamma(1,\frac1 \lambda)$
  - 可加性： $X_1,\cdots,X_n\sim Exp(\lambda)=\Gamma(1,\frac 1 \lambda),i.i.d.,\Sigma X_i\sim \Gamma(n,\frac 1 \lambda)$
  - 线性性： $X\sim \Gamma(1,\frac1 \lambda),nX\sim\Gamma(1,\frac1 {n\lambda})$
- $t$ 分布
  - $\xi\sim N(0,1),\eta\sim \chi^2(n)$ ，且二者独立，则 $T=\frac{\xi}{\sqrt{\eta/n}}\sim t(n)$ ，n个自由度的t分布


### 定义

- 无偏估计
  - 称 $\varphi(X_1 ,\cdots,X_n)$ 为 $g(\theta)$ 的无偏估计，若： $E_\theta\varphi(X_1,\cdots,X_n)=g(\theta) ,\forall\theta\in\Theta$
- 均方误差
  - $M_\theta(\varphi)=E_\theta[\varphi(X_1,\cdots,X_n)-g(\theta)]^2$
- 充分统计量
  - 称 $\varphi(X_1 ,\cdots ,X_n)$ 为 $g(\theta)$ 的充分统计量，若： $L(x_1,\cdots ,x_n;\theta)$ 可表示为 $q[\varphi(x_1 ,\cdots,x_n),g(\theta)]\cdot h(x)$
- 指数型分布
  - 称 $X$ 服从指数型分布，若 $f(x;\theta)=S(\theta)h(x)\exp\{\Sigma_{j=1}^k C_j(\theta)T_j(x)\}$
  - 此时似然函数为 $L(x_1,\cdots,x_n;\theta)=S(\theta)^n\prod h(x_i)\exp\{ \Sigma_{j=1}^k C_j(\theta)\Sigma_{i=1}^nT_j(x_i) \}$
  - 此时 $(\Sigma_{i=1}^nT_1(x_i),\cdots,\Sigma_{i=1}^nT_k(x_i))$ 是 $\theta$ 的一个充分统计量
- 完全统计量
  - 称 $\varphi(X_1 ,\cdots,X_n)$ 为 $g(\theta)$ 的完全统计量，若对任意函数 $u(\cdot)$ ，若 $E_\theta[u(\varphi(X_1,\cdots,X_n))]=0,\forall \theta\in\Theta$ ，则有 $u(\varphi(X_1,\cdots,X_n))\equiv0,\forall \theta\in\Theta$
- 一致最小方差无偏估计
  - 构造方法：完全充分统计量、C-R不等式
- C-R不等式
  - $Var_\theta(\varphi(X_1,\cdots,X_n))\geq\frac{(g(\theta)')^2}{nI(\theta)}$ ，其中Fisher信息量 $I(\theta)=E_\theta(\frac{\partial \log f(x;\theta)}{\partial\theta})^2$
  - 要求 $X$ 的支撑 $E=:\{ x:f(x;\theta)>0 \}$ 与 $\theta$ 无关


### 常用公式

- $\Sigma(x_i-\bar x^2) = \Sigma x_i^2 -n\bar x^2$
  

### 置信区间

- 已知 $\sigma^2=\sigma_0^2$ ，求 $\mu$ 置信区间
  - $X_i\sim N(\mu,\sigma_0^2),\bar X \sim N(\mu,\frac{\sigma_0^2}{n})$ 则 $\eta=:\frac{\bar X-\mu}{\sqrt{\sigma_0^2/n}}\sim N(0,1)$

    | $\gamma$ |$0.90$|$0.95$|$0.99$|
    | -| -| -| -|
    |$\eta$|$1.65$|$1.96$|$2.58$|

- 未知 $\sigma^2$ ，求 $\mu$ 置信区间
  - 用 $S^2$ 代替 $\sigma^2$ ， $T=:\frac{\bar X-\mu}{\sqrt{S^2/n}}$
  - 定理：
    - $X_1,\cdots,X_n\sim N(\mu,\sigma^2), i.i.d.,\frac{1}{\sigma^2}\Sigma_{i=1}^n(X_i-\bar X)^2\sim\chi^2(n-1)$
    - 上述 $\bar X$ 与 $\Sigma_{i=1}^n(X_i-\bar X)^2$ 独立
  - $T=\frac{\bar X-\mu}{\sqrt{S^2/n}}=\frac{\sqrt n(\bar X-\mu)/\sigma}{\sqrt{\frac{1}{(n-1)\sigma^2}\Sigma_{i=1}^n(X_i-\bar X)^2}}$ ，且分子分母独立，此时 $T=\frac{\bar X-\mu}{\sqrt{S^2/n}}\sim t(n-1)$
- 未知 $\mu$ ，求 $\sigma^2$ 置信区间
  - $Y=:\frac{(n-1)S^2}{\sigma^2}=\frac{\Sigma_{i=1}^n(X_i-\bar X)^2}{\sigma^2}\sim\chi^2(n-1)$
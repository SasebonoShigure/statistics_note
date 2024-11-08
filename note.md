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

## 第三章

### 定义

- 第I类错误： $H_0$ 真， $\varphi$ 否认 $H_0$
- 第II类错误： $H_0$ 假， $\varphi$ 不否认 $H_0$
- 功效函数 $\beta_\varphi(\theta)=:E_\theta\varphi(X)$
  - $\theta\in\Theta_0$ 时， 第I类错误概率为 $\beta_\varphi(\theta)$
  - $\theta\notin\Theta_0$ 时， 第II类错误概率为 $1-\beta_\varphi(\theta)$
- 称 $\varphi$ 是水平为 $\alpha$ 的检验法， 若 $\sup\{\beta_\varphi(\theta):\theta\in\Theta_0\}\leq\alpha$ ，即第I类错误不超过 $\alpha$ 
- 一致最大功效（UMP）检验：称 $\varphi$ 是水平为 $\alpha$ 的UMP检验，如果 $\varphi$ 水平为 $\alpha$ ，且对于任一检验水平为 $\alpha$ 的检验法 $\psi$ ， $\beta_\varphi(\theta)\geq\beta_\psi(\theta),\forall\theta\in\Theta\backslash\Theta_0$ ，即水平为 $\alpha$ 的检验法中第II类错误最小的
- 无偏检验：称 $\varphi$ 是水平为 $\alpha$ 的无偏检验，若对于 $\forall\theta\in\Theta\backslash\Theta_0$ ， $\beta_\varphi(\theta)\geq\alpha$
- 一致最大功效无偏（UMPU）检验： 无偏检验中第II类错误最小的
- 单参数指数型分布： 若 $X$ 的密度函数可表示为

  $$f(x,\theta)=S(\theta)h(x)e^{Q(\theta)V(x)}$$

  其中 $S(\theta)>0,h(x)>0,Q(\theta)$ 是 $\theta$ 的严格增函数，则称 $X$ 服从单参数指数型分布
  二项分布、指数分布、Poisson分布、一个参数已知的正态分布等皆为单参数指数型分布
  记 $t(X_1,\cdots,X_n)=t(X)=\Sigma_{i=1}^{n}V(X_i)$ 易知其为充分统计量

### 似然比检验法

- $\Theta=\{\theta_1,\theta_2\},\Theta_0=\{\theta_1\}$

  $$\begin{array}{lll}
    H_0:\theta=\theta_1 & \leftrightarrow & H_1:\theta=\theta_2
  \end{array}$$

  称 $\lambda(x)=L(x,\theta_2)/L(x,\theta_1)$ 为似然比
  - N-P引理：给定 $0\leq\alpha\leq 1$ ， 设

    $$W_0=\{x:L(x,\theta_2)>\lambda_0 L(x,\theta_1)\}$$

    $$\varphi_0(x)=\begin{cases}
      1 & \lambda(x)>\lambda_0(x) \\
      0 & 否则
    \end{cases}$$

    其中 $\lambda_0$ 满足 $E_{\theta_1}\varphi_0(x)=\beta_{\varphi_0}(\theta_1)=\alpha$ ， $\varphi_0$ UMP，无偏
- $X$ 服从单参数指数型分布
  - $\begin{array}{lll}
    H_0:\theta\leq\theta_1 & \leftrightarrow & H_1:\theta>\theta_2
    \end{array}$
    若存在 $C$ 满足 $P_{\theta_1}(\Sigma_{i=1}^{n}V(X_i))=\alpha$ 则UMP检验为

    $$\varphi_0(x)=\begin{cases}
      1 & \Sigma_{i=1}^{n}V(X_i)>C \\
      0 & 否则
    \end{cases}$$

  - $\begin{array}{lll}
    H_0:\theta\notin(\theta_1,\theta_2) & \leftrightarrow & H_1:\theta\in(\theta_1,\theta_2)
    \end{array}$
    若存在 $C_1<C_2$ 满足 $\begin{array}{ll}
    \beta_{\varphi_0}(\theta_i)=\alpha & (i=1,2)
    \end{array}$ 则UMP检验为

    $$\varphi_0(x)=\begin{cases}
      1 & C_1<\Sigma_{i=1}^{n}V(X_i)<C_2 \\
      0 & 否则
    \end{cases}$$

  - $\begin{array}{lll}
    H_0:\theta\in[\theta_1,\theta_2] & \leftrightarrow & H_1:\theta\notin[\theta_1,\theta_2]
    \end{array}$
    若存在 $C_1<C_2$ 满足 $\begin{array}{ll}
    \beta_{\varphi_0}(\theta_i)=\alpha & (i=1,2)
    \end{array}$ 则UMPU检验为

    $$\varphi_0(x)=\begin{cases}
      1 & \Sigma_{i=1}^{n}V(X_i)<C_1 或 \Sigma_{i=1}^{n}V(X_i)>C_2 \\
      0 & 否则
    \end{cases}$$

  - $\begin{array}{lll}
    H_0:\theta=\theta_0 & \leftrightarrow & H_1:\theta\neq\theta_0
    \end{array}$
    若存在 $C_1<C_2$ 满足
    1. $\beta_{\varphi_0}(\theta_0)=\alpha$
    2. $E_{\theta_0}(\varphi_0(X_1,\cdots,X_n)\Sigma_{i=1}^{n}V(X_i))=\alpha E_{\theta_0}(\Sigma_{i=1}^{n}V(X_i))$
    则UMPU检验为

    $$\varphi_0(x)=\begin{cases}
      1 & \Sigma_{i=1}^{n}V(X_i)<C_1 或 \Sigma_{i=1}^{n}V(X_i)>C_2 \\
      0 & 否则
    \end{cases}$$

    有对称性时，有时第一个条件蕴含第二个条件
  
### 广义似然比检验法

- 称 $\lambda(x)=\frac{\sup\{L(x,\theta):\theta\in\Theta\}}{\sup\{L(x,\theta):\theta\in\Theta_0\}}$ 为广义似然比。
- 广义似然比检验法：满足 $\sup\{\beta_{\varphi_0}(\theta):\theta\in\Theta_0\}=\alpha$ 的检验法

  $$\varphi_0(x)=\begin{cases}
    1 & \lambda(x)>\lambda_0(x) \\
    0 & 否则
  \end{cases}$$

  - 若 $t(x)$ 是
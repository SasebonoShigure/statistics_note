# 数理统计笔记

## 第二章

### 基本分布

- 两点（Bernoulli）分布， $X_i\sim B(1,p)$
  - $P(X_i=x_i)=p^{x_i}(1-p)^{1-x_i}$
- 指数分布， $X_i\sim Exp(\lambda)$
  - $p(x)=\lambda e^{-\lambda x}$
  - $E(x)=\frac{1}{\lambda},Var(x)= \frac{1}{\lambda ^2}$
  - MLE: $\hat\lambda = \frac{1}{\bar x}$
- 正态分布， $X\sim N(\mu,\sigma^2)$
  - $p(x) = \frac{1}{\sqrt{2\pi \sigma^2}}e^{-\frac{(x-\mu)^2}{2\sigma^2}}$
  - $E(x)=\mu, Var(x)=\sigma^2$
  - MLE: $\hat\mu=\bar x,\hat\sigma^2=\frac{1}{n}\sum(x_i-\bar x)^2$
  - 样本方差 $S^2=\frac{1}{n-1}\sum(x_i-\bar x)^2$ 是 $\sigma^2$ 的最小方差无偏估计
  - 性质：
    - $X\sim N(\mu,\sigma^2),nX\sim N(n\mu,n^2\sigma^2)$
    - $X_1,\cdots,X_n\sim N(\mu,\sigma^2), i.i.d.$

      $$\sum X_i\sim N(n\mu,n\sigma^2)$$

      $$\bar X \sim N(\mu,\frac{\sigma^2}{n})$$

- 均匀分布， $X\sim U[a,b]$
  - $E(x)=\frac{a+b}{2} , Var(x)=\frac{(b-a)^2}{12}$
  - MLE: $\hat a=\min\limits_{1\leq i\leq n} x_i,\hat b=\max\limits_{1\leq i\leq n} x_i$
  - 对于 $X_i \sim U[0,\theta]$ ，矩估计： $\tilde\theta=2\bar x$
- 泊松分布
  - $X\sim P(\lambda)=\frac{\lambda^k}{k!}e^{-\lambda}$
  - $E(X)=\lambda,Var(x)=\lambda$
  - $X_1,\cdots,X_n\sim P(\lambda),i.i.d.,\sum X_i\sim P(n\lambda)$
- 卡方分布
  - $X_1,\cdots,X_n\sim N(0,1), i.i.d.,\sum_{i=1}^nX_i^2\sim\chi^2(n)$
  - $\chi^2(n)=\Gamma(\frac n 2,\frac1 2)$
  - 性质：
    - $\xi\sim\chi^2(m), \eta\sim\chi^2(n)$ ，二者独立，

      $$\xi+\eta\sim\chi^2(m+n)$$

    - $X\sim\chi^2(n),EX=n$
- Gamma分布
  - 与指数分布关系： $X\sim Exp(\lambda)=\Gamma(1,\frac1 \lambda)$
  - 可加性： $X_1,\cdots,X_n\sim Exp(\lambda)=\Gamma(1,\frac 1 \lambda),i.i.d.$

    $$\sum X_i\sim \Gamma(n,\frac 1 \lambda)$$

  - 线性性： $X\sim \Gamma(1,\frac1 \lambda),$

    $$nX\sim\Gamma(1,\frac1 {n\lambda})$$

- $t$ 分布
  - $\xi\sim N(0,1),\eta\sim \chi^2(n)$ ，且二者独立，则

    $$T=\frac{\xi}{\sqrt{\eta/n}}\sim t(n)$$

    n个自由度的 $t$ 分布


### 定义

- 无偏估计
  - 称 $\varphi(X_1 ,\cdots,X_n)$ 为 $g(\theta)$ 的无偏估计，若： $E_\theta\varphi(X_1,\cdots,X_n)=g(\theta) ,\forall\theta\in\Theta$
- 均方误差
  - $M_\theta(\varphi)=E_\theta[\varphi(X_1,\cdots,X_n)-g(\theta)]^2$
- 充分统计量
  - 称 $\varphi(X_1 ,\cdots ,X_n)$ 为 $g(\theta)$ 的充分统计量，若： $L(x_1,\cdots ,x_n;\theta)$ 可表示为 $q[\varphi(x_1 ,\cdots,x_n),g(\theta)]\cdot h(x)$
- 指数型分布
  - 称 $X$ 服从指数型分布，若 $f(x;\theta)=S(\theta)h(x)\exp\{\sum_{j=1}^k C_j(\theta)T_j(x)\}$
  - 此时似然函数为 $L(x_1,\cdots,x_n;\theta)=S(\theta)^n\prod h(x_i)\exp\{ \sum_{j=1}^k C_j(\theta)\sum_{i=1}^nT_j(x_i) \}$
  - 此时 $(\sum_{i=1}^nT_1(x_i),\cdots,\sum_{i=1}^nT_k(x_i))$ 是 $\theta$ 的一个充分统计量，也是完全统计量。
- 完全统计量
  - 称 $\varphi(X_1 ,\cdots,X_n)$ 为 $g(\theta)$ 的完全统计量，若对任意函数 $u(\cdot)$ ，若 $E_\theta[u(\varphi(X_1,\cdots,X_n))]=0,\forall \theta\in\Theta$ ，则有 $u(\varphi(X_1,\cdots,X_n))\equiv0,\forall \theta\in\Theta$
- 一致最小方差无偏估计
  - 构造方法：完全充分统计量、C-R不等式
- C-R不等式
  - $Var_\theta(\varphi(X_1,\cdots,X_n))\geq\frac{(g(\theta)')^2}{nI(\theta)}$ ，其中Fisher信息量 $I(\theta)=E_\theta(\frac{\partial \log f(x;\theta)}{\partial\theta})^2$
  - 要求 $X$ 的支撑 $E=:\{ x:f(x;\theta)>0 \}$ 与 $\theta$ 无关


### 常用公式

- $\sum(x_i-\bar x)^2 = \sum x_i^2 -n\bar x^2$
- $\sum(x_i-\mu)^2 = \sum (x_i-\bar x)^2 +n(\bar x-\mu)^2$
- $X\sim N(0,1),E(X^4)=3;X\sim N(0,\sigma^2),E(X^4)=3\sigma^4$
  

### 置信区间

- 已知 $\sigma^2=\sigma_0^2$ ，求 $\mu$ 置信区间
  - $X_i\sim N(\mu,\sigma_0^2),\bar X \sim N(\mu,\frac{\sigma_0^2}{n})$ 则 $\eta=:\frac{\bar X-\mu}{\sqrt{\sigma_0^2/n}}\sim N(0,1)$

    | $\gamma$ |$0.90$|$0.95$|$0.99$|
    | -| -| -| -|
    |$\eta$|$1.65$|$1.96$|$2.58$|

- 未知 $\sigma^2$ ，求 $\mu$ 置信区间
  - 用 $S^2$ 代替 $\sigma^2$ ， $T=:\frac{\bar X-\mu}{\sqrt{S^2/n}}$
  - 定理：
    - $X_1,\cdots,X_n\sim N(\mu,\sigma^2), i.i.d.$

      $$\frac{1}{\sigma^2}\sum_{i=1}^n(X_i-\bar X)^2\sim\chi^2(n-1)$$

    - 上述 $\bar X$ 与 $\sum_{i=1}^n(X_i-\bar X)^2$ 独立
  - $T=\frac{\bar X-\mu}{\sqrt{S^2/n}}=\frac{\sqrt n(\bar X-\mu)/\sigma}{\sqrt{\frac{1}{(n-1)\sigma^2}\sum_{i=1}^n(X_i-\bar X)^2}}$ ，且分子分母独立，此时 

    $$T=\frac{\bar X-\mu}{\sqrt{S^2/n}}\sim t(n-1)$$

- 未知 $\mu$ ，求 $\sigma^2$ 置信区间

  $$Y=:\frac{(n-1)S^2}{\sigma^2}=\frac{\sum_{i=1}^n(X_i-\bar X)^2}{\sigma^2}\sim\chi^2(n-1)$$

## 第三章

### 定义

- 第I类错误
  - $H_0$ 真， $\varphi$ 否认 $H_0$
- 第II类错误
  - $H_0$ 假， $\varphi$ 不否认 $H_0$
- 功效函数 
  - $\beta_\varphi(\theta)=:E_\theta\varphi(X)$
    - $\theta\in\Theta_0$ 时， 第I类错误概率为 $\beta_\varphi(\theta)$
    - $\theta\notin\Theta_0$ 时， 第II类错误概率为 $1-\beta_\varphi(\theta)$
- 称 $\varphi$ 是水平为 $\alpha$ 的检验法， 若 $\sup\{\beta_\varphi(\theta):\theta\in\Theta_0\}\leq\alpha$ ，即第I类错误不超过 $\alpha$ 
- 一致最大功效（UMP）检验
  - 称 $\varphi$ 是水平为 $\alpha$ 的UMP检验，如果 $\varphi$ 水平为 $\alpha$ ，且对于任一检验水平为 $\alpha$ 的检验法 $\psi$ ， $\beta_\varphi(\theta)\geq\beta_\psi(\theta),\forall\theta\in\Theta\backslash\Theta_0$ ，即水平为 $\alpha$ 的检验法中第II类错误最小的
- 无偏检验
  - 称 $\varphi$ 是水平为 $\alpha$ 的无偏检验，若对于 $\forall\theta\in\Theta\backslash\Theta_0$ ， $\beta_\varphi(\theta)\geq\alpha$
- 一致最大功效无偏（UMPU）检验
  - 无偏检验中第II类错误最小的
- 单参数指数型分布
  - 若 $X$ 的密度函数可表示为

    $$f(x,\theta)=S(\theta)h(x)e^{Q(\theta)V(x)}$$

    其中 $S(\theta)>0,h(x)>0,Q(\theta)$ 是 $\theta$ 的严格增函数，则称 $X$ 服从单参数指数型分布
    记 $t(X_1,\cdots,X_n)=t(X)=\sum_{i=1}^{n}V(X_i)$ 易知其为充分统计量
  - 二项分布、指数分布、Poisson分布、一个参数已知的正态分布等皆为单参数指数型分布
- F分布
  - $\xi\sim\chi^2(n_1),\eta\sim\chi^2(n_2)$ ，且相互独立，则 $\zeta=\frac{\xi/n_1}{\eta/n_2}\sim F(n_1,n_2)$
  - 性质
    - $X\sim F(n_1,n_2)$ ，则

      $$\frac 1 X \sim F(n_2,n_1)$$

    - $T\sim t(n)$ ，则

      $$T^2 \sim F(1,n)$$

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
    若存在 $C$ 满足 $P_{\theta_1}(\sum_{i=1}^{n}V(X_i))=\alpha$ 则UMP检验为

    $$\varphi_0(x)=\begin{cases}
      1 & \sum_{i=1}^{n}V(X_i)>C \\
      0 & 否则
    \end{cases}$$ 

  - $\begin{array}{lll}
    H_0:\theta\notin(\theta_1,\theta_2) & \leftrightarrow & H_1:\theta\in(\theta_1,\theta_2)
    \end{array}$
    若存在 $C_1<C_2$ 满足 $\begin{array}{ll}
    \beta_{\varphi_0}(\theta_i)=\alpha & (i=1,2)
    \end{array}$ 则UMP检验为

    $$\varphi_0(x)=\begin{cases}
      1 & C_1<\sum_{i=1}^{n}V(X_i)<C_2 \\
      0 & 否则
    \end{cases}$$

  - $\begin{array}{lll}
    H_0:\theta\in[\theta_1,\theta_2] & \leftrightarrow & H_1:\theta\notin[\theta_1,\theta_2]
    \end{array}$
    若存在 $C_1<C_2$ 满足 $\begin{array}{ll}
    \beta_{\varphi_0}(\theta_i)=\alpha & (i=1,2)
    \end{array}$ 则UMPU检验为

    $$\varphi_0(x)=\begin{cases}
      1 & \sum_{i=1}^{n}V(X_i)<C_1 或 \sum_{i=1}^{n}V(X_i)>C_2 \\
      0 & 否则
    \end{cases}$$

  - $\begin{array}{lll}
    H_0:\theta=\theta_0 & \leftrightarrow & H_1:\theta\neq\theta_0
    \end{array}$
    若存在 $C_1<C_2$ 满足
    1. $\beta_{\varphi_0}(\theta_0)=\alpha$
    2. $E_{\theta_0}(\varphi_0(X_1,\cdots,X_n)\sum_{i=1}^{n}V(X_i))=\alpha E_{\theta_0}(\sum_{i=1}^{n}V(X_i))$
    则UMPU检验为

    $$\varphi_0(x)=\begin{cases}
      1 & \sum_{i=1}^{n}V(X_i)<C_1 或 \sum_{i=1}^{n}V(X_i)>C_2 \\
      0 & 否则
    \end{cases}$$

    有对称性时，有时第一个条件蕴含第二个条件
  
### 广义似然比检验法

- 广义似然比
  - $$\lambda(x)=\frac{\sup\{L(x,\theta):\theta\in\Theta\}}{\sup\{L(x,\theta):\theta\in\Theta_0\}}$$
- 广义似然比检验法
  - 满足 $\sup\{\beta_{\varphi_0}(\theta):\theta\in\Theta_0\}=\alpha$ 的检验法

    $$\varphi_0(x)=\begin{cases}
      1 & \lambda(x)>\lambda_0(x) \\
      0 & 否则
    \end{cases}$$

  - 若 $t(x)$ 是 $\theta$ 的充分统计量，则 $\lambda(x)$ 可表示为 $t(x)$ 的函数
- $X\sim N(\mu,\sigma^2)$ ， $\sigma^2$ 未知， $\begin{array}{lll}
    H_0:\mu=\mu_0 & \leftrightarrow & H_1:\mu\neq\mu_0
  \end{array}$
  - $\lambda(x)$ 表示为 $|T|$ 的严格单调增函数， $\mu=\mu_0$ 时 $T=\frac{\bar X-\mu}{\sqrt{S^2/n}}\sim t(n-1)$

    $$\varphi_0(x)=\begin{cases}
      1 & |T|>C \\
      0 & 否则
    \end{cases}$$

    $P_{\mu_0\sigma}(|T|>C)=\alpha$

- $X\sim N(\mu,\sigma^2)$ ， $\sigma^2$ 未知， $\begin{array}{lll}
    H_0:\mu\leq\mu_0 & \leftrightarrow & H_1:\mu>\mu_0
  \end{array}$
  - 同上
- $X\sim N(\mu,\sigma^2)$ ， $\mu$ 未知， $\begin{array}{lll}
    H_0:\sigma^2=\sigma_0^2 & \leftrightarrow & H_1:\sigma^2\neq\sigma_0^2
  \end{array}$
  - 令 $u = \frac{1}{\sigma^2}\sum_{i=1}^n(x_i-\bar x)^2$ ， $\lambda(x)=f(u)=e^{-\frac n 2}n^{\frac n 2}u^{-\frac n 2}e^{\frac u 2}$

    $$\varphi_0(x)=\begin{cases}
      1 & u<C_1 或 u>C_2 \\
      0 & 否则
    \end{cases}$$

    由 $u$ 分布得出
- $X\sim N(\mu,\sigma^2)$ ， $\mu$ 未知， $\begin{array}{lll}
    H_0:\sigma^2\leq\sigma_0^2 & \leftrightarrow & H_1:\sigma^2>\sigma_0^2
  \end{array}$

  - $$\varphi_0(x)=\begin{cases}
      1 & u>C_2 \\
      0 & 否则
    \end{cases}$$

- 两个正态总体 $X\sim N(\mu_1,\sigma_1^2),Y\sim N(\mu_2,\sigma_2^2)$ ， $X,Y$ 独立，样本为 $x=(x_1,\cdots,x_{n_1}),y=(y_1,\cdots,y_{n_2})$

  1. $$\begin{array}{lll}
      H_0:\sigma_1^2=\sigma_2^2 & \leftrightarrow & H_1:\sigma_1^2\neq\sigma_2^2
      \end{array}$$

      - 令 $F=\frac{S_1^2}{S_2^2}$ ， 是 $S_1^2,S_2^2$ 无偏估计之比，且当 $H_0$ 成立时 $F\sim F(n_1-1,n_2-1)$

      $$\varphi_0(x,y)=\begin{cases}
        1 & F<C_1 或 F>C_2 \\
        0 & 否则
      \end{cases}$$

  2. $$\begin{array}{lll}
      H_0:\sigma_1^2\leq\sigma_2^2 & \leftrightarrow & H_1:\sigma_1^2>\sigma_2^2
      \end{array}$$

      $$\varphi_0(x,y)=\begin{cases}
        1 & F>C_2 \\
        0 & 否则
      \end{cases}$$

  3. $\sigma_1^2,\sigma_2^2$ 未知，已知 $\sigma_1^2=\sigma_2^2$ ，

      $$\begin{array}{lll}
        H_0:\mu_1=\mu_2 & \leftrightarrow & H_1:\mu_1\neq\mu_2
      \end{array}$$

      - $$T=\frac{\bar x-\bar y}{\sqrt{\sum(x_i-\bar x)^2+\sum(y_j-\bar y)^2}}\sqrt{\frac{n_1n_2(n_1+n_2-2)}{n_1+n_2}}$$

        $\lambda(x,y)$ 可表示为 $|T|$ 的单调增函数。当 $H_0$ 成立时 $T\sim t(n_1+n_2-2)$

        $$\varphi_0(x,y)=\begin{cases}
          1 & |T|>C \\
          0 & 否则
        \end{cases}$$

        $C$ 满足 $\mu_1=\mu_2$ 时 $P(|T|>C)=\alpha$
  4. $\sigma_1^2,\sigma_2^2$ 未知，已知 $\sigma_1^2=\sigma_2^2$ ，

      $$\begin{array}{lll}
        H_0:\mu_1\leq\mu_2 & \leftrightarrow & H_1:\mu_1>\mu_2
      \end{array}$$

      - $$\varphi_0(x,y)=\begin{cases}
          1 & T>C_2 \\
          0 & 否则
        \end{cases}$$

        $C$ 满足 $\mu_1=\mu_2$ 时 $P(|T|>C_2)=\alpha$

### 拟合优度检验

- $\chi^2$ 检验
  - $H_0$ ：$X$ 来自某个类型的分布
  - 先MLE估计出 $\hat\theta$ ，再在此条件下做检验，检验问题变为 $\tilde H_0 :\theta=\hat\theta$
  - 离散情况，共 $m+1$ 个取值，连续情况将数轴分为 $m+1$ 段
  - $v_k$ ：落入第k段样本频数， $p_k$ ：理论落入第k段概率

    $$V=:\sum_{k=1}^{m+1}(\frac{v_k}{n}-p_k)^2\frac n {p_k}=\sum_{k=1}^{m+1}\frac{(v_k-np_k)^2}{np_k}$$

    近似地， $V\sim \chi^2(m)$

## 第四章

### 一元线性回归

- 模型为 $Y=a+bx+\varepsilon$
  - 数据为

    $$\begin{cases}
    y_1=a+bx_1+\varepsilon_1 \\
    \cdots \\
    y_n=a+bx_n+\varepsilon_n
    \end{cases} $$

  - 对模型的假设由弱到强

    $$\begin{array}{lrcl}
        1. & E(\varepsilon_i)=0 &\leftrightarrow& E(Y_i)=ax_i+b \\
        2. & Var(\varepsilon_i)=\sigma^2 &\leftrightarrow& Var(Y_i)=\sigma^2 \\
        3. & \varepsilon_i 相互独立，&\leftrightarrow& Y_i 相互独立，\\
           & Cov(\varepsilon_i,\varepsilon_j)=0(i\neq j)& & Cov(Y_i,Y_j)=0(i\neq j)\\
        4. & \varepsilon_i\sim N(0,\sigma^2) &\leftrightarrow& Y_i\sim N(a+bx_i,\sigma^2)
      \end{array}$$

- 最小二乘估计
  - 使

    $$Q(a,b)=\sum_{i=1}^{n}(y_i-a-bx_i)^2=\sum_{i=1}^{n}\varepsilon_i$$

    达到最小的 $\hat{a},\hat{b}$ 称为 $a,b$ 的LSE
    求偏导，解得

    $$ \begin{cases}
        \hat{a}=\bar{y}-\hat{b}\bar{x} &(直线过(\bar x,\bar y))\\
        \hat{b}=\frac{\sum(x_i-\bar x)(y_i-\bar y)}{\sum(x_i-\bar x)^2}
      \end{cases} $$

    记 $l_{xx}=\sum(x_i-\bar x)^2,l_{xy}=\sum(x_i-\bar x)(y_i-\bar y)$  ，则 $b=l_{xy}/l_{xx}$
- 线性关系检验：
  - 平方和分解公式
    - 残差： $e_i=y_i-\hat y_i$
      偏差平方和： $l_{yy}=\sum (y_i-\bar y)^2$
      回归平方和： $U=\sum (\hat y_i-\bar y)^2$
      残差平方和： $Q=\sum (y_i-\hat y_i)^2=\sum e_i^2$
      $Q=Q(\hat a,\hat b)$
  - 引理1： $l_{yy}=U+Q$


  - $$\begin{array}{rcl}
      H_0:b=0 & \leftrightarrow & H_1:b\neq 0
      \end{array} $$

    - 若误差 $\varepsilon_i\sim N(0,\sigma^2)$ 且 $H_0:b=0$ 成立时，

      $$F=:\frac{U}{Q/(n-2)}\left( =\frac{\frac{U}{\sigma^2}/1}{\frac{Q}{\sigma^2}/(n-2)} \right)\sim F(1,n-2)$$

      给定$\alpha$ 查表 $F>\lambda$ 时否认 $H_0$

  - $\hat{a},\hat{b}$ 是 $a,b$ 的无偏估计。
    当 $\varepsilon_i\sim N(0,\sigma^2)$ 成立时 $\hat{a},\hat{b}$ 是MLE，
    - 此时 $\sigma^2$ 的MLE为 $Q(\hat a,\hat b)/n$ ， 无偏估计为 $Q/(n-2)$
  - 相关系数

    $$r=\frac{l_{xy}}{\sqrt{l_{xx}l_{yy}}} $$

    - $|r|$ 越接近 $1$ ， $xy$ 线性相关程度越强
    - $r^2=\frac{U}{U+Q}$ ，即回归平方和占偏差平方和的比例
    - $F=\frac{U}{Q/(n-2)}=\frac{r^2}{(1-r^2)/(n-2)}$ ，关于 $r^2$ 单增
- 预测
  - $Y_0$ 点估计 $\hat Y_0=\hat a+\hat b x_0$ ，也是 $E(Y_0)$ 的点估计
  - $Y_0$ 置信区间：

    $$T=\frac{Y_0-\hat Y_0}{\sqrt{dQ/(n-2)}}= \frac{Y_0-\hat a-\hat b x_0}{\sqrt{dQ/(n-2)}}\sim t(n-2)$$

    其中 $d=1+\frac 1 n + \frac{(x_0-\bar x)^2}{l_{xx}}$

- 控制
  - 给定 $1-\alpha$ ，希望 $Y_0$ 在 $[A,B]$ 内
    - $x=x_0$ 时， $Y_0$ 的 $1-\alpha$ 水平置信区间为

      $$[\hat Y_0-\lambda\sqrt{dQ/(n-2)},\hat Y_0+\lambda\sqrt{dQ/(n-2)}]$$

      让此置信区间始终在 $[A,B]$ 内，找满足条件的 $x$

### 线性模型

- 模型为 $Y=\beta_1x_1+\beta_2x_2+\cdots\beta_px_p+\varepsilon$
  - 数据为

    $$\begin{cases}
    y_1=\beta_1x_{11}+\cdots+\beta_px_{1p} \\
    \cdots \\
    y_n=\beta_1x_{n1}+\cdots+\beta_px_{np}
    \end{cases} $$

  - 记

    $$Y = \begin{pmatrix}
      Y_1 \\
      \vdots \\
      Y_n
      \end{pmatrix}_{n \times 1}, \quad
      X = \begin{pmatrix}
      x_{11} & \cdots & x_{1p} \\
      \vdots & \ddots & \vdots \\
      x_{n1} & \cdots & x_{np}
      \end{pmatrix}_{n \times p}, \quad
      \beta = \begin{pmatrix}
      \beta_1 \\
      \vdots \\
      \beta_p
      \end{pmatrix}_{p \times 1}, \quad
      e = \begin{pmatrix}
      e_1 \\
      \vdots \\
      e_n
      \end{pmatrix}_{n \times 1}$$

    模型可表示为 $Y=X\beta+e$
  - 对 $e$ 有两种假设：

    $\begin{array}{lll}
        A: & E(e)=0, & Cov(e,e)=E(e,e')=\sigma^2I
      \end{array}$

    假设零均值，不相关，同方差；不假设独立

    $\begin{array}{ll}
        B: & e\sim N(0,\sigma^2I)
      \end{array}$

    假设正态，且 $e_i$ 互相独立(正态不相关等价于独立)

- 最小二乘估计
  - 记$Q(\beta) = \sum_{i=1}^n e_i^2 = e' \cdot e = \sum_{i=1}^n \left[ y_i - (\beta_1 x_{i1} + \cdots + \beta_p x_{ip}) \right]^2= (Y - X\beta)'(Y - X\beta) = \| Y - X\beta \|^2$
    最小化 $Q(\beta)$ 的 $\hat\beta$ 为LSE
  - 定理3.1
      1. 满足 $X'X\beta=X'Y$ 的 $\hat\beta$ 是LSE(正规方程总有解但不一定唯一)，一定存在
      2. LSE唯一的充要条件是 $X$ 满秩，此时 $\hat\beta=(X'X)^{-1}X'Y$
  - 定理3.2 设 $X$ 满秩，且假定A成立
      1. $E(\hat\beta)=\beta$
      2. $Cov(\hat{\beta}, \hat{\beta}) = \sigma^2 (X'X)^{-1}$
      3. $EQ(\hat{\beta}) = (n - p)\sigma^2$ (不满秩时为$n-r,r=rank(X)$)， $\sigma^2$ 无偏估计为 $Q(\hat{\beta})/(n-r)$

- 线性可估性
  - 设 $c=(c_1,\cdots,c_p)'\in R^p$ 为任意p维向量
    定义：称 $c'\beta$ 是（线性）可估的，若存在 $a=(a_1,\cdots ,a_n)'\in R^n$ ，使得 $E(a'Y)=c'\beta$ 对任意的 $\beta\in R^p$ 成立
    - $X$ 满秩时任意 $c\in R^p$ ， $c'\beta$ 都可估，此时 $a'=c'(X'X)^{-1}X'$ ； $X$ 不满秩时只有部分 $c'\beta$ 可估
  - 定理3.3
    - 设假定A成立，则 $c'\beta$ 可估的充要条件为 $c\in\mu(X')$ ，即 $c'$ 可表示为 $X$ 的行向量的线性组合。
  - 定理3.4
    - 设假定A成立， $c'\beta$ 线性可估，则任取 $\hat\beta$ 为 $\beta$ 的一个LSE， $c'\hat\beta$ 必为 $c'\beta$ 的唯一的最小方差线性无偏估计。
  - 定义：若 $\hat\beta$ 是 $\beta$ 的LSE，则称 $c'\hat\beta$ 为 $c'\beta$ 的LSE。
  - 定理3.5
    - 设假定A成立， $X$ 的秩为 $r$ ，则 $EQ(\hat{\beta}) = (n-r)\sigma^2$ ，即 $\hat\sigma^2=Q(\hat{\beta})/(n-r)$ 是 $\sigma^2$ 的无偏估计。
    - 假定B成立时， $\hat\sigma^2$ 是 $\sigma^2$ 的一致最小方差无偏估计
- 带约束的线性模型
  - 约束：$H\beta=r_0$ ，其中 $H:s\times p$ ，秩为s
  - $\hat\beta$ 是满足约束且使 $\|Y-X\beta\|^2$ 达到最小的 $\beta$
  1. 消去多余参数
  2. 拉格朗日乘子
      - 设未知量 $\lambda\in R^s$ ，解方程

        $$\begin{cases}
          X'X\hat{\beta} - H'\lambda = X'Y \\
          H\hat{\beta} = r_0
          \end{cases}$$

- 参数线性相关性检验
  - $H_0 : H\beta = 0 \leftrightarrow H_1 : H\beta \neq 0$ ， $H: s\times p$ ，每一行是一个待检验的假设，假定B成立
    - 记 $W=\mu(x),$ 维数为 $r,W_0=\{\eta=X\beta ,H\beta = 0\}$ ，秩为 $0<q<r$ ，
      $\Theta=\{(\beta,\sigma^2)\},\Theta_0=\{(\beta,\sigma^2),H\beta=0\}$ ， $\hat\beta,\hat\beta_0$ 分别为无约束和有约束的LSE
    - $\hat{\xi} = Proj_{W}(Y) = X\hat{\beta} ,\hat{\xi}_0 = Proj_{W_0}(Y) = X\hat{\beta}_0,\| Y - \hat{\xi} \|^2=Q$
    - 广义似然比检验法， $L(Y, \beta, \sigma^2) = (2\pi\sigma^2)^{-\frac{n}{2}} e^{-\frac{1}{2\sigma^2} \| Y - X\beta \|^2}$ 

      $$F = \frac{\| \hat{\xi} - \hat{\xi}_0 \|^2 / (r - q)}{\| Y - \hat{\xi} \|^2 / (n - r)} = \frac{\| X\hat{\beta} - X\hat{\beta}_0 \|^2 / (r - q)}{\| Y - X\hat{\beta} \|^2 / (n - r)}$$

      广义似然比 $\lambda$ 是 $F$ 的单调增函数，在假定B成立前提下，若 $H_0$ 成立， $F\sim F(r-q,n-r)$
      检验法：

      $$\varphi_0=\begin{cases}
          1 & F>C \\
          0 & 否则
        \end{cases}$$

    - 在检验哪些 $\beta$ 可以删去时，$r-q$ 变为删去的 $\beta$ 数量 $k$ ($k=rank(H)$)
    - 此处有正交分解 $\| Y - \hat{\xi}_0 \|^2=\| Y - \hat{\xi} \|^2+\| \hat{\xi} - \hat{\xi}_0 \|^2$ ， $\| \hat{\xi} - \hat{\xi}_0 \|^2=\tilde Q -Q$
    - 推论： $X\hat\beta$ 与 $\hat e=Y-X\hat\beta,Q$ 独立， $Q/\sigma^2\sim \chi^2(n-r)$
- 参数线性组合的检验与置信区间
  - 称 $a$ 为 $c$ 的伴随元，若 $a'Y=c'\hat\beta$ 是 $c'\beta$ 无偏估计。 $X$ 满秩时 $a=X(X'X)^{-1}c$
  - LSE 的 $c'\hat\beta$ 与 $Q$ 独立
  - 设 $c'\beta$ 可估，记 $\hat\sigma^2=Q/(n-r)$

    $$\frac{c'(\hat{\beta} - \beta)}{\hat{\sigma} \| a \|} \sim t(n - r)$$

    其中 $c'\beta=E(Y_0)$ ，(可以求$E(Y_0)置信区间$)
    $X$ 满秩时即：

    $$\frac{c'(\hat{\beta} - \beta)}{\hat{\sigma} \sqrt{c'(X'X)^{-1}c}} \sim t(n - p)$$
  - 对于假设 $H_0 : c'\beta = r_0$ ，检验统计量为：

    $$\frac{c'\hat{\beta}-r_0}{\hat{\sigma} \| a \|} \sim t(n - r)$$

  - 预测 $Y_0$ ，统计量为

    $$T = \frac{c'\hat{\beta} - Y_0}{\hat{\sigma} \sqrt{\| a \|^2 + 1}} \sim t(n - r)$$

    预测时 $c=x_0$

## 第五章

### 全面试验的方差分析

- 单因素试验的方差分析
  - 单因素 $A$ 可取 $A_1,\cdots,A_s$ 个水平，每个水平重复 $r$ 次。
    模型为 $Y_{ij}=\mu_i+e_{ij},i=\{1,\cdots,s\},j=\{1,\cdots,r\}$
    假设 $e_{ij}:i.i.d.,e_{ij}\sim N(0,\sigma^2)$ ， $\sigma^2$ 未知，共 $s+1$ 个未知参数
  - 检验问题为 $H_0:\mu_1=\cdots=\mu_s$
  - $S_T=S_e+S_A$
    - $S_T = \sum_{i=1}^s \sum_{j=1}^r \left(Y_{ij} - \bar{Y}\right)^2$
    - $S_e=\sum_{i=1}^s \sum_{j=1}^r \left(Y_{ij} - \bar{Y}_{i\cdot}\right)^2$
    - $S_A=r \sum_{i=1}^s \left(\bar{Y}_{i\cdot} - \bar{Y}\right)^2$
  - $H_0$ 成立时 $S_A$ 相对于 $S_e$ 较小，检验统计量为

    $$F = \frac{S_A/(s - 1)}{S_e/s(r - 1)}\sim F(s-1,s(r-1))$$
- 两因素试验的方差分析
  - 因素 $A$ 可取 $A_1,\cdots,A_s$ 个水平，$B$ 可取 $B_1,\cdots,B_t$ 个水平，每个水平重复 $r$ 次。
  - 模型为 $Y_{ijk}=\mu_{ij}+e_{ijk},i=\{1,\cdots,s\},j=\{1,\cdots,t\},k=\{1,\cdots,r\}$
  - 假设 $e_{ijk}:i.i.d.,e_{ijk}\sim N(0,\sigma^2)$ ， $\sigma^2$ 未知，共 $st+1$ 个未知参数
  - $S_T=S_e+S_A+S_B+S_{A\times B}$
    - $S_e=\sum_{ijk} \left( Y_{ijk} - \overline{Y}_{ij.} \right)^2$
    - $S_{A\times B}=r \sum_{ij} \left( \overline{Y}_{ij.} - \overline{Y}_{i..} - \overline{Y}_{.j.} + \overline{Y} \right)^2$
    - $S_A=t \cdot r \sum_{i} \left( \overline{Y}_{i..} - \overline{Y} \right)^2 $
    - $S_B=s \cdot r \sum_{j} \left( \overline{Y}_{.j.} - \overline{Y} \right)^2$
    - $S_A,S_B,S_{A\times B}$ 分别为 $A,B$ 主效应和 $AB$ 的交互作用
  - 统计量为

    $$F_1 = \frac{S_A/(s-1)}{S_e/st(r-1)} \sim F(s-1, st(r-1))$$

    $$F_2 = \frac{S_B/(t-1)}{S_e/st(r-1)} \sim F(t-1, st(r-1))$$

    $$F_3 = \frac{S_{A \times B}/(s-1)(t-1)}{S_e/st(r-1)} \sim F((s-1)(t-1), st(r-1))$$

## 第七章

### 统计决策问题

- 定义
  - 决策函数： $\delta(X_1,\cdots,X_n)$
  - 损失函数： $L(\theta,a)$ 表示参数真值为 $\theta$ 时采取决策 $a$ 的损失
  - $R(\theta,\delta)=E_\theta L(\theta,\delta(X_1,\cdots,X_n))$ 为决策 $\delta$ 的风险函数
  - 最优决策 $\delta^*$ ： $R(\theta, \delta^*) \leq R(\theta, \delta) \quad \forall \theta \in \Theta$
  - $\delta$ 是可容许的，若不存在另一 $\delta'$ 使 $R(\theta, \delta') \leq R(\theta, \delta) \quad \forall \theta \in \Theta$ 且对某一 $\theta_0$ 严格成立
  - $minimax$ 决策： $\sup\{R(\theta, \delta^*), \theta \in \Theta\} \leq \sup\{R(\theta, \delta), \theta \in \Theta\}$

### 贝叶斯统计

- 基本公式
  - $f(x, \theta) = L(x \mid \theta) \pi(\theta) = \pi(\theta \mid x) h(x)$
    - 后验分布 $\pi(\theta \mid x)$ 正比于 $L(x \mid \theta) \pi(\theta)$
    - $h(x)=\int_{\Theta} L(x \mid \theta) \pi(\theta) \, d\theta$
  - 贝叶斯置信区间 $[\theta_L,\theta_U]$
    - $\int_{\theta_L}^{\theta_U} \pi(\theta \mid x) \, d\theta = 1 - \alpha$
    - 可去两边各为 $\alpha/2$
  - 预测
    - $f(x_{n+1} \mid x_1, \cdots, x_n) = \int_{\Theta} f(x_{n+1} \mid \theta) \cdot \pi(\theta \mid x_1, \cdots, x_n) \, d\theta$
  - $\delta$ 的平均风险： $\rho(\delta) = E_\pi R(\theta, \delta) = \int_{\Theta} R(\theta, \delta) \pi(\theta) \, d\theta$
    - 使 $\rho(\delta)$ 达到最小的 $\delta^*$ 为贝叶斯决策
    - 实际中找最优决策：找行动 $a_x$ 使 $\int_{\Theta} L(\theta, \delta(x)) \pi(\theta \mid x) \, d\theta$ 达到最小， $\delta^*(x)=a_x$
- 共轭分布族
  - 定义：如果先验分布 $\pi(\theta)$ 的类型使得后验分布 $\pi(\theta \mid \mathbf{x})$ 仍为此类型，则称 $\pi(\theta)$ 是 $f(\mathbf{x} \mid \theta)$ 的共轭分布。
  - 设 \(X_1, \cdots, X_n\) 是来自正态分布 \(N(\mu_0, \sigma^2)\) 的简单随机样本，其中 \(\mu_0\) 已知，若取参数 \(\sigma^2\) 的先验分布为 \(IG(\alpha, \beta)\)，则 \(\sigma^2\) 的后验分布为：
    \[
    IG\left(\alpha + \frac{n}{2}, \, \beta + \frac{1}{2} \sum_{i=1}^n (X_i - \mu_0)^2\right)
    \]

    - 后验期望是第二个参数除以第一个参数

  
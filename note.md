# 数理统计笔记

## 第二章

### 基本分布

- 两点（Bernoulli）分布， $X_i\sim B(n,p)$
  - $P(X_i=x_i)=p^{x_i}(1-p)^{1-x_i}$
- 指数分布，$Xi\sim Exp(\lambda)$
  - $p(x)=\lambda e^{-\lambda x}$
  - $ E(x)=\frac{1}{\lambda},Var(x)= \frac{1}{\lambda ^2}$
  - MLE: $ \hat\lambda = \frac{1}{\bar x} $
- 正态分布，$ X\sim N(\mu,\sigma^2) $
  - $ p(x) = \frac{1}{\sqrt{2\pi \sigma^2}}e^{-\frac{(x-\mu)^2}{2\sigma^2}} $
  - $ E(x)=\mu, Var(x)=\sigma^2 $
  - MLE: $ \hat\mu=\bar x,\hat\sigma^2=\frac{1}{n}\Sigma(x_i-\bar x)^2 $
  - 样本方差 $S^2=\frac{1}{n-1}\Sigma(x_i-\bar x)^2$ 是 $\sigma^2$ 的最小方差无偏估计
- 均匀分布，$ X\sim U[a,b] $
  - $ E(x)=\frac{1}{b-a} , Var(x)=\frac{(b-a)^2}{12} $
  - MLE: $ \hat a=\min\limits_{_1\leq i\leq n} x_i,\hat b=\max\limits_{_1\leq i\leq n} x_i $
  - 对于$ X_i \sim U[0,\theta] $，矩估计：$ \tilde\theta=2\bar x $

### 定义

- 无偏估计
  - 称 $\phi(X_1 ,\cdots,X_n) $ 为 $ g(\theta) $ 的无偏估计，若：$ E_\theta\phi(X_1,\cdots,X_n)=g(\theta) ,\forall\theta\in\Theta $
- 均方误差
  - $ M_\theta(\phi)=E_\theta[\phi(X_1,\cdots,X_n)-g(\theta)]^2 $
- 充分统计量
  - 称 $ \phi(X_1 ,\cdots,X_n) $ 为 $ g(\theta) $ 的充分统计量，若：$L(x_1,\cdots,x_n;\theta)$ 可表示为 $q[\phi(x_1 ,\cdots,x_n),g(\theta)]\cdot h(x)$
- 指数型分布
  - 称 $X$ 服从指数型分布，若 $ f(x;\theta)=S(\theta)h(x)\exp\{ \Sigma_{j=1}^k C_j(\theta)T_j(x) \} $
  - 此时似然函数为$ L(x_1,\cdots,x_n;\theta)=S(\theta)^n\prod h(x_i)\exp\{ \Sigma_{j=1}^k C_j(\theta)\Sigma_{i=1}^nT_j(x_i) \} $
  - 此时 $ (\Sigma_{i=1}^nT_1(x_i),\cdots,\Sigma_{i=1}^nT_k(x_i)) $ 是 $\theta$ 的一个充分统计量
- 完全统计量
  - 称 $ \phi(X_1 ,\cdots,X_n) $ 为 $ g(\theta) $ 的完全统计量，若对任意函数 $u(\cdot)$ ，若 $ E_\theta[u(\phi(X_1,\cdots,X_n))]=0,\forall \theta\in\Theta $，则有 $u(\phi(X_1,\cdots,X_n))\equiv0,\forall \theta\in\Theta$
- 一致最小方差无偏估计
  - 构造方法：完全充分统计量、C-R不等式
- C-R不等式
  - $ Var_\theta(\phi(X_1,\cdots,X_n))\geq\frac{(g(\theta)')^2}{nI(\theta)} $，其中Fisher信息量 $ I(\theta)=E_\theta(\frac{\partial \log f(x;\theta)}{\partial\theta})^2 $
  - 要求 $X$ 的支撑 $ E=:\{ x:f(x;\theta)>0 \} $ 与 $\theta$ 无关
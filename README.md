This project attempts to implement the adaptive Lasso using a method prescribed by Rajaratnam et al (2016). This method relies on the foundation set by Park & Casella (2008) and uses their conditionals as a starting point.
Conditional probability distributions are used to express uncertainty in the model coefficients for variables of interest. Rajaratnam et al found that along with expressing uncertainty, the structure of the conditionals can also be exploited to estimate the Lasso point estimate itself. 

For computing this estimate, they suggested shrinking the Gibbs sampler to the limit of σ2 -> 0. Also, at this limit, the conditional sequence of the Gibbs sampler reduces down to a deterministic sequence. This formalized technique is called Shrinkage via the Limit Of Gibbs sampler, or SLOG.

The modification adds a 10-fold Leave One Out Cross Validation component to the original algorithm.


# Brent's recollections about AI use

Claude suggested Lemma lem:specialization should be written like this:  (I've only made minor changes)

  \begin{lemma}
  \label{lem:specialization}
  Let $(S^=, S^{\ne})$ be a differentially simple system computed with a
  constant-lowest ranking, and let $c^* \in \mathbb{C}^n$ lie in the system's
  cell.  Then
  \begin{enumerate}
  \item[\rm(i)] each $p \in S^=$ and its specialization $p(c^*)$ have the same
    leader, to the same degree;
  \item[\rm(ii)] $H_{S^=(c^*)} = \{\, g(c^*) : g \in H_{S^=} \,\}$, and no element
    of it is zero; and
  \item[\rm(iii)] $(S^=(c^*), S^{\ne}(c^*))$ is again differentially simple.
  \end{enumerate}
  \end{lemma}

  Because $c^*$ lies in the cell, \cite{ThomasDecomp}, Remark~2.3 applied upward
  from the constants supplies a solution $\sigma$ of $(S^=, S^{\ne})$ whose
  constant coordinates are $c^*$.  Let $p \in S^=$ have leader $u$.  By
  \condref{A2} the initial $I_p$ does not vanish at $\sigma$, so $I_p(c^*)$ is not
  the zero polynomial --- were it identically zero it would vanish at $\sigma$ in
  particular.  Hence $\deg_u p(c^*) = \deg_u p$, the leader of $p(c^*)$ is still
  $u$, and its initial is $I_p(c^*)$: this is (i) and half of (ii).  For the
  separant, $\partial / \partial u$ acts on the coefficients of the expansion of
  $p$ in $u$ while the substitution $c \mapsto c^*$ touches no occurrence of $u$,
  so the two commute,
  \[ 
  s_{p(c^*)} \;=\; \partial\, p(c^*) / \partial u \;=\; (\partial p / \partial u)(c^*) \;=\; s_p(c^*) ,
  \]
  and this too is non-zero, not vanishing at $\sigma$ by
  Lemma~\ref{lem:initials-and-separants-specialize}.


but I've made some changes: it didn't define H_{S^=} in the statement of the lemma, so I added that.

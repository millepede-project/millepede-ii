\file
Draft Manual by V.Blobel.

No code - only documentation for doxygen.

\page draftman_page Millepede II - Draft Manual

**Linear Least Squares Fits with a Large Number of Parameters**

\author V. Blobel, University Hamburg, 2007

\remark
Adapated to formatting with Doxygen (C. Kleinwort, DESY, 2012). The abstract
has been moved to section \ref intro_sec. Some clarifications as result of user
feedback included (C. Kleinwort, DESY, 2013). A list of major changes to the
implementation described in this manual can be found \ref changes_page "here".
The solution method GMRES described later has never been implemented.
Instead \ref ch-minres "MINRES" can be used.

\tableofcontents

\section ssec_preface Preface

<b>Linear least squares problems.</b>
The most important general method for the determination of
parameters from measured data is the linear least squares method.
It is usually stable and accurate, requires no initial values
of parameters and is
able to take into account all the correlations between different
parameters, even if there are many parameters.

<b>Global and local parameters.</b>
The parameters of a least squares problem can sometimes be distinguished
as *global* and *local* parameters. The mathematical model underlying
the measurement depends on both types of parameters. The interest is
in the determination of the global parameters, which appear in all
the measurements. There are many sets of local parameters, where
each set appears only in one, sometimes small, subset of measurements.

<b>Alignment of large detectors in particle physics.</b>
Alignment problems for large detectors in particle physics often require the
determination of a large number of alignment parameters, typically of the
order of \f$100\f$ to \f$1\,000\f$, but sometimes above \f$10\,000\f$.
Alignment parameters for example define the
accurate space coordinates and orientation of detector components.
In the alignment usually special alignment measurements are combined with
data of particle reaction, typically tracks from physics interactions and
from cosmics. In this paper the alignment parameters are called *global*
parameters. Parameters of a single
track like track slopes and curvatures are called *local* parameters.

One approximate alignment method is to perform least squares fits on the data
e.g. of  single tracks assuming fixed alignment parameters.
The *residuals*, the *deviations* between
the fitted and measured data, are then used to estimate the alignment
parameters afterwards. Single fits depend only on the
small number of local parameters
like slope or curvature of the track, and are easy to solve.
This approximate method however is not a correct method, because the
(local) fits depend on a wrong model, they ignore the global
parameters, and the result are biased fits. The adjustment of parameters
based on the observed (biased) residuals will then result in
biased alignment parameters. If these alignment parameters are
applied as corrections in repeated fits, the remaining
residuals will be reduced, as desired. However, the fitted parameters
will still be biased. In practice this procedure is often
applied iteratively; it is however not clear wether the procedure
is converging.

A more efficient and faster method is an overall least squares fit,
with all the global
parameters and local parameters, perhaps from thousands or millions of events,
determined simultaneously. The **Millepede** algorithm makes use
of the special structure of the least squares matrices in such a
simultaneous fit: the  global parameters can be determined from a
matrix equation, where the matrix dimension is given by the
number of global parameters only, irrespective of the total number of
local parameters, and without any approximation.
If \f$n\f$ is the number of global parameters, then an
equation with a symmetric \f$n\f$-by-\f$n\f$ matrix has to be solved. The
**Millepede** programm has been used for up to \f$n \approx 5000\f$
parameters in the past. Solution of the matrix equation was done
with a fast program for the solution of matrix equations with
inversion of the symmetric matrix, which on  a standard PC takes
a time \f$t \approx 2 \times  10^{-8} \, \textrm{sec}\, \times n^3\f$. Even for
\f$n = 5000\f$ the computing time is, with
\f$t = 2 \times 10^{-8} \, \textrm{sec}\, \times
5000^3 = 40\f$ minutes, still acceptable.

The next-generation detectors in particle physics however require up to
about 100000 global parameters and with the previous solution method the
space- and time-consumption is too large
for a standard PC. Memory space is needed for the full symmetric
matrix, corresponding to \f$1/2 n^2 \times 8\f$ bytes for double precision, if the
symmetry of the matrix is observed, and if solution is done *in-space*.
This means memory space of 400 Mbyte and
40 Gbyte for \f$n = 10\,000\f$ and \f$n = 100\,000\f$, and computing times of
about 6 hours and almost a year, respectively.

<b>Millepede I and II.</b>
The second version of the **Millepede** programm, described here,
allows to use several different methods for the solution of the matrix equation
for global parameters, while keeping the algorithm, which decouples the
determination of the local parameters, i.e. still taking
into account all correlations between all global parameters.
The matrix can be stored as a sparse matrix, if
a large fraction of off-diagonal elements has
zero content.
A complete matrix inversion of a large sparse matrix would
result in a full matrix, and would take an unacceptable time. Different
methods can be used in **Millepede II** for a sparse matrix.
The symmetric matrix can be
accumulated in a special sparse storage scheme.
Special solution methods can be used, which require only products
of the symmetric sparse matrix with different vectors.
Depending on the fraction of
vanishing matrix elements, solutions for up to a number of 100000
global parameters should be possible on a standard PC.

The structure of **Millepede II** is different from the **Millepede I**
version. The accumulation of data and the solution are splitted,
with the solution in a stand-alone program. Furthermore the
global parameters are now characterized by a label (any positive
integer) instead of a (continuous) index. Those features should
simplify the application of the program for alignment problems
even for large number of global parameters in the range of \f$10^5\f$.

\section sec_math Mathematical Methods

\subsection ssec_large Large problems with global and local parameters

The solution  of optimization problems requires the determination
of certain parameters. In certain optimization problems the parameters
can be classified as either *global* or *local*.
This classification can be used, when the input data of the
optimization problem
consist out of a potentially large group of data sets, where the
description of each single set of data requires certain *local*
parameters, which appear only in the description of one data set.
The *global* parameter may appear in all data sets, and the
interest is in the determination of these global parameters. The
local parameters can be called nuisance parameters.

An example for an optimization problem with global and local parameters
from experimental high energy physics is the alignment of track detectors
using a large number of track data. The track data are position
measurements of a charged particle track, which can be fitted for example
by a helix with five parameters, if the track detector is in a homogeneous
magnetic field. The position measurement, from several detector planes
of the track detector, depend on the position and orientation of certain
sensors, and the corresponding coordinates and angles are global
parameters. In an alignment the values of the global parameters are
improved by minimizing the deviations between measurement and
parametrization of a large number of tracks. Numbers of tracks
in the order of one Million, with of the order of 10 data poinst per track,
and of \f$1000\f$ to \f$100,000\f$ global parameters are typical for
modern track detectors in high energy physics.

Optimal values of the global parameters require a simultaneous fit
of all parameters with all data sets. A straight-forward ansatz for such a
fit seems to require to fit Millions of parameters, which is impossible.
An alternative is to perform each local fit separately, fitting the
local parameters only. From the residuals of the local fits one can then
try to fit the values of the global parameters. This approach neglects
the correlations between the global and local parameters. Nevertheless
the method can be applied iteratively with the hope, that the
convergence is not to slow and does not require too many iterations.

The optimization problem is complicated by the fact, that often not all
global parameters are defined by the ansatz. In the alignment example
the degrees of freedom describing a translation and rotation of the whole
detector are not fixed. The solution of the optimization problem requires
therefore certain  equality constraints, for example zero overall
translation and rotation of the detector; these constraints are described
by a linear combination of global parameters. Such an equality
constraint can only be satisfied in an overall fit, not with separated
local and global fits.
\note Constraints in optimization
are either equality or inequality constraints, which have exactly
to be taken into account in the solution. The term *constraint*
is often used in a loose sense, like: ``the parameters are constrained
by our measurement''.

Global parameters are denoted by the vector \f$\vec{p}\f$ and local parameters
for the data set \f$j\f$ by the vector \f$\vec{q}_j\f$. The objective function
to be minimized in the global parameter optimization is the sum of
the local objective Function, depending on the global parameters \f$\vec{p}\f$
and the vectors \f$\vec{q}_j\f$:
\f{equation*}{
        F(\vec{p}, \vec{q}) = \sum_j F_j(\vec{p}, \vec{q}_j)
\f}
In the least squares method the objective function to be minimized is
a sum of the squares of residuals \f$z_i\f$ between the measured values and the
parametrization, weighted by the inverse variance of the measured value:
\f{equation*}{
        F_j(\vec{p}, \vec{q}_j) = \tfrac{1}{2}
\sum_i \frac{z_i^2}{\sigma_i^2}
\f}
with the residual \f$z_i = y_i - f_i(\vec{p},\vec{q}_j)\f$, where
\f$y_i\f$ is the measured value and \f$f_i(\vec{p},\vec{q}_j)\f$ is the
corresponding parametrization.
This form of the objective function assumes independent measurements,
with a single variance value \f$\sigma_i\f$ assigned.
\note
In experimental high energy physics the objective function is usually called a
\f$\chi^2\f$-function. In statistics there is a \f$\chi^2\f$-distribution with a
well-defined meaning. The
minimum value \f$\times 2\f$  of the objective function follows,
under certain conditions,
the \f$\chi^2\f$-distribution.

\subsection ssec_opt Optimization without constraints

\subsubsection ssec_qm The quadratic model

The standard optimization method for smooth objective functions is the
**Newton** method. The objective function \f$F(\vec{p})\f$, depending on
a parameter vector \f$\vec{p}\f$, is approximated by a quadratic model
\f$\tilde{F}(\vec{p})\f$
\f{equation*}{ \label{eq:qapp}
\tilde{F}_k \left( \vec{p}_k + \vec{d} \right)
= F_k + \vec{g}^{\top}  \vec{d} + \tfrac{1}{2}
\vec{d}^{\top} \mathbf{C}  \vec{d}
\f}
where \f$\vec{p}_k\f$ is the vector \f$\vec{p}\f$ in the \f$k\f$-th iteration and
where the vector \f$\vec{g}\f$ is the gradient of the objective function;
the matrix \f$\mathbf{C}\f$ is the Hessian (second derivative matrix) of the
objective function \f$F(\vec{p})\f$ or an approximation of the Hessian.
The minimum of the quadratic approximation requires the gradient
to be equal to zero.
A step \f$\vec{d}\f$ in parameter space is calculated by the solution of the
matrix equation, obtained from the derivative of the quadratic model:

\anchor eq-cdmg (1)
\f{equation*}{  \label{eq:cdmg}
                \mathbf{C} \, \vec{d} = - \vec{g} \, .
\f}
With the correction vector \f$\vec{d}\f$ the new value \f$\vec{p}_{k+1} =
\vec{p}_k + \vec{d}\f$ for the next iteration is obtained. The
matrix \f$\mathbf{C}\f$ is a constant in a *linear least squares* problem
and the minimum is determined in a single step (no iteration necessary).

\subsubsection sssec_par Partitioning of matrices

The special structure of the matrix \f$\mathbf{C}\f$ in a matrix equation
\f$\mathbf{C} \vec{d} = - \vec{g}\f$
may allow a significant simplification of the solution.
Below the symmetric matrix \f$\mathbf{C}\f$  is partitioned
into submatrices, and
the vectors \f$\vec{d}\f$ and \f$\vec{g}\f$  are partitioned into two subvectors;
then the matrix equation can be written in the form
\f{equation*}{
\left(
    \begin{array}{ccc|c}
&           & &         \\
&    \mathbf{C}_{11} & & \mathbf{C}_{21}^{\top}  \\
&           & &         \\   \hline
&     \mathbf{C}_{21} & & \mathbf{C}_{22}
    \end{array}
    \right)
\left( \begin{array}{c}
~\\   \vec{d_1}  \\  ~\\  \hline \vec{d_2}
        \end{array} \right)  = -
\left( \begin{array}{c}
~\\   \vec{g_1}  \\  ~\\  \hline  \vec{g_2}
        \end{array} \right)  \; ,
\f}
where the submatrix \f$\mathbf{C}_{11}\f$ is a
\f$p\f$-by-\f$p\f$ square
matrix and the submatrix \f$\mathbf{C}_{22}\f$ is a \f$q\f$-by-\f$q\f$
square matrix, with \f$p+q=n\f$, and \f$\mathbf{C}_{21}\f$ is a \f$q\f$-by-\f$p\f$ matrix.
Now it is assumed that the inverse of the \f$q\f$-by-\f$q\f$
sub-matrix \f$\mathbf{C}_{22}\f$ is available. In certain problems this
may be easily calculated, for example if \f$\mathbf{C}_{22}\f$ is diagonal.

If the sub-vector \f$\vec{d}_1\f$ would not exist, the solution for the
sub-vector \f$\vec{d}_2\f$ would be defined by the matrix equation
\f$\mathbf{C}_{22} \; \vec{d}_2^{*} = - \vec{g}_2\f$,
where the star indicates the special character of this solution, which is

\anchor eq-spa2 (2)
\f{equation*}{ \label{eq:spa2}
            \vec{d_2}^{*} = -  \mathbf{C}_{22}^{-1} \;  \vec{g_2}  \; .
\f}

Now, having the inverse sub-matrix \f$\mathbf{C}_{22}^{-1}\f$, the submatrix
of the complete inverse matrix \f$\mathbf{C}\f$ correponding to the upper left
part \f$\mathbf{C}_{11}\f$ is the inverse of the
symmetric \f$p\f$-by-\f$p\f$ matrix

\anchor eq-Schur (3)
\f{equation*}{ \label{eq:Schur}
\mathbf{S} =   \mathbf{C}_{11} - \mathbf{C}_{21}^{\top} \mathbf{C}_{22}^{-1}
    \mathbf{C}_{21}     \; ,
\f}
the so-called *Schur complement*.
With this matrix \f$\mathbf{S}\f$ the solution of the whole matrix equation can
be written in the form
\f{equation*}{  \label{eq:solvea1}
\left( \begin{array}{c}
~\\   \vec{d_1}  \\  ~\\  \hline \vec{d_2}
        \end{array} \right) = -
\left(
    \begin{array}{ccc|c}
&           & &         \\
&    \mathbf{S}^{-1} &  & - \mathbf{S}^{-1} \mathbf{C}_{21}^{\top}
            \mathbf{C}_{22}^{-1}  \\
&           & &         \\   \hline
&   - \mathbf{C}_{22}^{-1} \mathbf{C}_{21} \mathbf{S}^{-1}
& & \mathbf{C}_{22}^{-1} - \mathbf{C}_{22}^{-1}
\mathbf{C}_{21} \mathbf{S}^{-1} \mathbf{C}_{21}^{\top} \mathbf{C}_{22}^{-1}
    \end{array}
    \right)
\left( \begin{array}{c}
~\\   \vec{g_1}  \\  ~\\  \hline  \vec{g_2}
        \end{array} \right)   \; .
\f}

The sub-vector \f$\vec{d_1}\f$ can be obtained from the solution of the
matrix equation

\anchor eq-solvea2 (4)
\f{equation*}{ \label{eq:solvea2}
\left(
    \begin{array}{ccc}
& & \\ &  \mathbf{C}_{11} - \mathbf{C}_{21}^{\top} \mathbf{C}_{22}^{-1}
    \mathbf{C}_{21} \\
& & \end{array} \right)
\left( \begin{array}{c}
~\\   \vec{d_1}  \\  ~\\
        \end{array} \right) = -
\left( \begin{array}{c}
~\\   \vec{g_1} -  \mathbf{C}_{21}^{\top} \mathbf{C}_{22}^{-1} \vec{g}_2   \\  ~\\
        \end{array} \right)
=-
\left( \begin{array}{c}
~\\   \vec{g_1} -  \mathbf{C}_{21}^{\top}  \vec{d_2}^*   \\
~\\
        \end{array} \right)
\f}
using the known right-hand-side of the equation with
the special solution \f$\vec{d_2}^*\f$.

In a similar way the
vector \f$\vec{d_2}\f$ could be calculated. However,
if the interest is the determination of this sub-vector \f$\vec{d_1}\f$
only, while
the sub-vector \f$\vec{d_2}\f$ is not needed, then only the  equation
\ref eq-solvea2 "(4)" has to be solved after calculation of
the special solution \f$\vec{d_2}^*\f$ (equation \ref eq-spa2 "(2)") and the Schur
complement \f$\mathbf{S}\f$ (equation \ref eq-Schur "(3)").
Some computer time can be saved by this
method, especially if the matrix \f$\mathbf{C}_{22}^{-1}\f$ is easily calculated
or already known before; note that the matrix \f$\mathbf{C}_{22}\f$ does not
appear directly in the solution, only the inverse
\f$\mathbf{C}_{22}^{-1}\f$.

This method  of removing unnecessary parameters
was already known in the nineteenth century.
The method can be applied repeatedly,
and therefore may simplify the solution of problems with a
large number of parameters.
The method is not an approximation, but is exact and it takes into account
all the correlations introduced by the removed parameters.
\note One example is: Schreiber, O. (1877): Rechnungsvorschriften f&uuml;r die
trigonometrische Abteilung der Landesaufnahme, Ausgleichung und
Berechnung der Triangulation zweiter Ordnung. Handwritten notes. Mentioned
in W. Jordan (1910): Handbuch der Vermessungskunde, Sechste erw. Auflage,
Band I, Paragraph III: 429-433. J.B.Metzler, Stuttgart.

\subsubsection sssec_locpar Local parameters

A set of local measured data \f$y_i\f$ is considered.
The local data \f$y_i\f$ are assumed to be
described by a linear or non-linear
function \f$f(x_i,\vec{q})\f$,
depending on a (small) number of local parameters \f$\vec{q}\f$.
\f{equation*}{ \label{eq:measur}
        y_i = f(x_i,\vec{q}) + \varepsilon_i  \; .
\f}
The parameters \f$\vec{q}\f$ are called the  local parameters,
valid for the specific group of measurements (local-fit object).
The quantity \f$\varepsilon\f$ is the measurement error, with standard deviation
\f$\sigma_i\f$.
The quantity \f$x_i\f$ is assumed to be the coordinate of the measured value
\f$y_i\f$, and is one argument of the function \f$f(x_i,\vec{q})\f$.

The local parameters are  determined in a least squares fit.
If the function
\f$f(x_i,\vec{q})\f$ depends *non-linearly* on the local parameters
\f$\vec{q}\f$, an iterative procedure is used, where the function
is linearized, i.e. the *first* derivatives of the function
\f$f(x_i,\vec{q})\f$ with respect to the local parameters \f$\vec{q}\f$
are calculated. The function is thus expressed as a linear function of local
parameter corrections \f$\vec{\Delta q}\f$ at some reference value \f$\vec{q}_k\f$:

\anchor eq-ydeff (5)
\f{equation*}{ \label{eq:ydeff}
        f(x_i,\vec{q}_k +\vec{\Delta q})  =
        f(x_i,\vec{q}_k) +
        \frac{\partial f}{\partial q_1} \Delta q_1
    + \frac{\partial f}{\partial q_2} \Delta q_2 + \ldots  \; ,
\f}
where the derivatives are calculated for \f$ \vec{q} \equiv \vec{q}_k\f$.
For each single measured value,
the residual measurement \f$z_i\f$
\f{equation*}{
z_i \equiv  y_i - f(x_i,\vec{q}_k)
\f}
is calculated.
For each iteration a linear system of equations (normal equations of least
squares) has to be solved for the parameter corrections \f$\vec{\Delta q}\f$
with a matrix \f$\vec{\Gamma}\f$ and a gradient
vector \f$\vec{g}\f$ with elements

\anchor eq-normaloc (6)
\f{equation*}{ \label{eq:normaloc}
\Gamma_{jk} = \sum_i \left( \frac{\partial f_i}{\partial q_j} \right)
        \left( \frac{\partial f_i}{\partial q_k} \right)
            \frac{1}{\sigma_i^2}
            \quad \quad \quad \quad
\beta_j= \sum_i \left( \frac{\partial f_i}{\partial q_j} \right)
\,  \frac{z_i}{\sigma_i^2}
\; ,
\f}
where the sum is over all measurements \f$y_i\f$ of the local-fit object.
Corrections \f$\vec{\Delta q}\f$ are determined by the solution of the
matrix equation
\f{equation*}{
            \vec{\Gamma} \vec{\Delta q} = - \vec{\beta} \; ,
\f}
and a new reference value is obtained by
\f{equation*}{
        \vec{q}_{k+1} = \vec{q}_k + \vec{\Delta q}
\f}
and then, with \f$k\f$ increased by 1, this is repeated until convergence
is reached.

\subsubsection sssec_glopar Global parameters

Now global parameters are considered, which contribute to all the
measurements. The expectation function \ref eq-ydeff "(5)"
is extended to include corrections for
global parameters. Usually only few of the global parameters
influence a local-fit object.
A global parameter is identified by a label \f$\ell\f$;
assuming that labels \f$\ell\f$ from a set
\f$\Omega\f$ contribute to a single measurement the extended equation
becomes

\f{equation*}{ \label{eq:zdefey}
        z_i = y_i - f(x_i,\vec{q},\vec{p}) =
\sum_{j=1}^{\nu} \left( \frac{\partial f}{\partial q_j} \right)  \Delta q_j
+ \sum_{\ell \in \Omega}
\left( \frac{\partial f}{\partial p_{\ell}} \right) \Delta p_{\ell}  \; .
\f}


\subsubsection sssec_simfit The simultaneous fit of global and local parameters

In the following it is assumed that there is a set of \f$N\f$ local measurements.
Each local measurement, with index \f$j\f$, depends on \f$\nu\f$ local parameters
\f$\vec{q}_j\f$, and all of them depend on the global parameters.
In a simultaneous fit of all global parameters plus local parameters
from \f$N\f$ subsets of the data there are in
total \f$(n+N\cdot\nu)\f$ parameters, and the standard solution requires the
solution of \f$(n+N\cdot\nu)\f$ equations with a computation proportional to
\f$(n+N \cdot \nu)^3\f$. In the next chapter it is shown, that the
problem can be reduced to a system of \f$n\f$ equations, for the global
parameters only.

For a set of  \f$N\f$ local measurements one obtains a system of
least squares normal equations with large
dimensions, as is shown in equation \ref eq-huge "(7)"
The matrix on the left side of equation \ref eq-huge "(7)"
has, from each local measurement, three types of contributions.
The first part is a contribution of a symmetric matrix \f$\vec{C_1}_j\f$, of
dimension \f$n\f$ (number of global parameters), and is calculated from the
(global) derivatives \f$\partial f/\partial p_{\ell}\f$.
All the matrices  \f$\vec{C_1}_j\f$
are added up in the upper left corner of the big matrix of the normal
equations. The second contribution is the symmetric matrix
\f$\vec{\Gamma}_j\f$ (compare equation \ref eq-normaloc "(6)"),
which gives a contribution to the big matrix on the
diagonal and is depending only on the \f$j\f$-th local measurement and the
(local) derivatives \f$\partial f/\partial q_j\f$.
The third (mixed) contribution is a rectangular matrix \f$\mathbf{G}_j\f$, with
a row number of \f$n\f$ (global) and a column number of \f$\nu\f$ (local).
There are two contributions to the vector of the normal equations (gradient),
\f$\vec{g_1}_j\f$ for the global and \f$\vec{\beta}_j\f$ for the local parameters.
The complete matrix equation is given by

\anchor eq-huge (7)
\f{equation*}{  \label{eq:huge}  \renewcommand{\arraystretch}{1.2}
\left(
    \begin{array}{ccc||ccc|c|ccc}
&           & &  & & & & & &         \\
&   \sum \vec{C_1}_j  & &  & \cdots & &  \mathbf{G}_j & & \cdots &  \\
&           & &  & & & & & &        \\   \hline \hline
&   & &  &  & & & & & \\
&  \vdots & & & \ddots & & 0 & & 0 & \\
&   & & & &  & & & & \\  \hline
&    \mathbf{G}^{\top}_j &  &  & 0 & & \vec{\Gamma}_j & &0 &  \\ \hline
&   & & & & & &  & & \\
&  \vdots & & & 0 & & 0 & & \ddots & \\
&   & & & & & & & &  \\
    \end{array}
    \right)
. \left( \begin{array}{c}
\\  \vec{d} \\  \\  \hline  \hline
\\ \vdots \\  \\  \hline
\vec{\Delta q}_j \\ \hline
\\ \vdots \\  \\
        \end{array} \right)
=
- \left( \begin{array}{c}
\\  \sum \vec{g_1}_j \\  \\  \hline \hline
\\ \vdots \\  \\  \hline
\vec{\beta}_j \\ \hline
\\ \vdots \\  \\
        \end{array} \right)
\f}
In this matrix equation the matrices  \f$\vec{C_1}_j\f$,
\f$\vec{\Gamma}_j\f$, \f$\mathbf{G}_j\f$
and the vectors  \f$\vec{g_1}_j\f$ and \f$\vec{\beta}_j\f$
contain contributions from the \f$j\f$-th local measurement.
Ignoring the global parameters (i.e. keeping them constant)
one could solve the normal equations
\f$ \vec{\Gamma}_j \vec{\Delta q}_j^* = - \vec{\beta}_j\f$
for each local measurement separately by
\f{equation*}{  \label{eq:ignore}
\vec{\Delta q}_j^* = - \vec{\Gamma}_j^{-1} \vec{\beta}_j \, .
\f}
The complete system of normal equations has a special structure, with many
vanishing sub-matrices. The only connection between the local parameters of
different partial measurements is given by the sub-matrices
\f$\mathbf{G}_j\f$ und \f$\vec{C_1}_j\f$,

\subsubsection sssec_redsize Reduction of matrix size

The aim of the fit is solely to determine the global parameters;
final best parameters of the local parameters are not needed.
The matrix of equation \ref eq-huge "(7)" is written in a partitioned form.
The general solution can also be written in partitioned form.
Many of the sub-matrices of the huge matrix in equation \ref eq-huge "(7)"
are zero and this has the effect, that the
formulas for the sub-matrices of the inverse matrix are very simple.

By this procedure the \f$n\f$ normal equations

\anchor eq-nsb (8)
\f{equation*}{   \label{eq:nsb}   \renewcommand{\arraystretch}{1.2}
\left(
    \begin{array}{ccc}
&           &      \\
&    \mathbf{C}  &  \\
&           &      \\
    \end{array}
    \right)
\left( \begin{array}{c}
\\  \vec{d} \\  \\
        \end{array} \right)
= -
\left( \begin{array}{c}
\\  \vec{g}  \\  \\
        \end{array} \right)   \; ,
\f}
are obtained, which only contain the global parameters, with a
modified matrix \f$\mathbf{C}\f$ and a modified vector \f$\vec{g}\f$,
\f{equation*}{  \label{eq:nsc}
\mathbf{C} =  \sum_j \vec{C_1}_j + \sum_j \vec{C_2}_j
\quad \quad  \quad \quad
\vec{g} =    \sum_j \vec{g_1}_j + \sum_j  \vec{g_2}_j
\f}
with the following local contributions to  \f$\mathbf{C}\f$ and  \f$\vec{g}\f$
from the \f$j\f$-th local fit:

\anchor eq-nsc2 (9)
\f{equation*}{  \label{eq:nsc2}
\vec{C_2}_j =  - \mathbf{G}_j \vec{\Gamma}_j^{-1} \mathbf{G}_j^{\top}
\quad \quad \quad \quad
\vec{g_2}_j =
- \mathbf{G}_j \left( \vec{\Gamma}_j^{-1} \vec{\beta}_j\right)
= - \mathbf{G}_j \vec{\Delta q}_j^*   \; .
\f}
The set of normal equations \ref eq-nsb "(8)" contains explicitly only the global
parameters; implicitly it contains, through the correction matrices,
the complete information from the local parameters, influencing the
fit of the global parameters. The parentheses in equation \ref eq-nsc2 "(9)"
represents the solution for the local parameters, ignoring the global
parameters.
The solution
\f{equation*}{ \label{eq:ared}
\vec{d} = - \mathbf{C}^{-1}\,  \vec{g}
\f}
represents the solution vector \f$\vec{d}\f$ with covariance matrix
\f$\mathbf{C}^{-1}\f$.
The solution is direct, no iterations or approximations are required.
The dimension of the matrix to compute \f$\vec{d}\f$ from equation
\ref eq-cdmg "(1)" is reduced from  \f$(n+N\cdot\nu)\f$ to \f$n\f$. The vector
\f$\vec{d}\f$ is the correction for the global parameter vector \f$\vec{p}\f$.
Iterations may be necessary for other reasons, namely
* the equations depend *non-linearly* on the global
parameters; the equations have to be linearized;
* the data contain outlier, which have to be removed in a sequence
of cuts, becoming narrower during the iteration, or which have to
be down-weighted;
* the accuracy of the data is not known before, and has to be
determined from the data (after the alignment).

For iterations the vector
\f$\vec{d}\f$ is the correction for the global parameter vector \f$\vec{p}_k\f$
in iteration \f$k\f$ to obtain the global parameter vector \f$\vec{p}_{k+1}
= \vec{p}_k + \vec{d}\f$ for the next iteration.

\subsubsection sssec_nonlin Nonlinear least squares

A method for *linear* least squares fits with a large number
of parameters, perhaps with *linear* constraints, is discussed
in this paper. Sometime of course the model is *nonlinear*
and also constraints may be nonlinear. The standard method to treat
these problems is linearization: the nonlinear equation is replaced
by a linear equation for the correction of a parameter (Taylor
expansion); this requires a good approximate value of the parameter.
In principle this method requires an iterative improvement of the
parameters, but sometimes even one iteration may be sufficient.

The treatment of nonlinear equations is not directly supported by
the program package, but it will in general not be too difficult to
organize a program application with nonlinear equations.

\subsubsection sssec-outlow Outliers

Cases with very large residuals within the data can distort the result
of the method of least squares. In the method of M-estimates, the
least squares method is modified to Maximum likelihood method.
Basis is the residual between
data and function value (expected data value), normalized by the
standard deviation:
\f{equation*}{
        \zeta = \frac{y - f(x)}{\sigma}
\f}
In the method of M-estimates the objective function \f$F(.)\f$
to be minimized is defined in terms of a probability density function
of the normalized residual
\f{equation*}{
        F(.)  = \sum_i \rho(\zeta_i)
            \quad \quad \quad \quad
    \rho(\zeta) = \ln \textrm{pdf}(\zeta)
\f}
From the probability density function \f$\textrm{pdf}(\zeta)\f$
a *influence function* \f$\psi(\zeta)\f$
is defined
\f{equation*}{
\textrm{influence function} \; \psi(\zeta) = \text{d} \rho(\zeta)/\text{d}\zeta
        \quad \quad \quad \quad
\textrm{additional weight factor} \; \omega(\zeta) =  \psi(\zeta)/\zeta
\f}
and a weight in addition to the normal least squares weight
\f$w_i = 1/\sigma_i^2\f$ can be derived from the influence function for
the calculation of the (modified) normal equations of least squares.

For the standard least squares method the function \f$\rho(\zeta)\f$ is simply
\f$\rho(\zeta) =   \zeta^2/2\f$
and it follows, that the influence function is \f$\psi(\zeta)=\zeta\f$
and the weight
factor is \f$\omega(\zeta) = 1\f$ (i.e. no extra weight).
The influence function value
increases with the value of the normalized residual without limits, and thus
outliers have a large and unlimited influence.
In order to reduce the influence of outliers, the probability density
function \f$\textrm{pdf}(\zeta)\f$ has to be modified for large values of
\f$|\zeta|\f$ to avoid the unlimited increase of the influence. Several
functions \f$\textrm{pdf}(\zeta)\f$ are proposed,
for example the Huber function and
the Cauchy function.

<b> Huber function:</b> A simple function is the
Huber function, which is quadratic and thus identical to least squares for
small \f$|\zeta|\f$,
but linear for larger  \f$|\zeta|\f$, where the influence function becomes
a constant \f$C_{\textrm{H}} \cdot \textrm{sign}(\zeta)\f$:
\f{equation*}{ \label{eq:huber}
\textrm{Huber function: pdf} \quad
\rho(\zeta) = \begin{cases} \zeta^2/2 \\
C_{\textrm{H}} \left( |\zeta| - C_{\textrm{H}}/2 \right) \\
    \end{cases}
\quad  \quad
\textrm{factor} \quad
\omega(\zeta) = \begin{cases} 1
&\textrm{if} \; |\zeta| \le  C_{\textrm{H}}  \\
C_{\textrm{H}}/|\zeta|   &\textrm{if} \;
|\zeta| >  C_{\textrm{H}} \end{cases}
\f}
The extra weight is
\f$\omega(\zeta) = C_{\textrm{H}}/|\zeta|\f$ for large  \f$|\zeta|\f$.
A standard value for
\f$C_{\textrm{H}}\f$ is \f$C_{\textrm{H}}=1.345\f$;
for this value the efficiency for Gaussian data without outliers
is still 95 \%.
For very large deviations the additional weight factor decreases with
\f$1/|\zeta|\f$.

<b> Cauchy function:</b>
For small deviation the Cauchy function is close to the least squares
expression, but for large deviations it increases only logarithmically.
For very large deviations the additional weight factor decreases with
\f$1/\zeta^2\f$:
\f{equation*}{ \label{eq:cauchy}
\textrm{Cauchy function: pdf} \quad
\rho(\zeta) =  \tfrac{1}{2} C_{\textrm{c}} \ln \left( 1 +
\left( \zeta/ C_{\textrm{c}}\right)^2 \right)
\quad \quad \quad
\textrm{factor} \quad
\omega = 1/\left( 1 + \left( \zeta/ C_{\textrm{c}} \right)^2 \right)
\f}
A standard value is \f$C_{\textrm{c}}=2.3849\f$;
for this value the efficiency for Gaussian data without outliers
is still 95 \%.

\subsection ssec-lincon Optimization with linear constraints

The minimization of an objective function \f$F(\vec{p})\f$ is
often not sufficient.
Several degrees of freedom may be undefined and require
additional conditions, which can be expressed as equality constraints.
For \f$m\f$ linear equality constraints the problem is the
minimization of a non-linear function \f$F(\vec{p})\f$ subject to a set of
linear constraints:
\f{equation*}{  \label{eq:cproblem}
\min F(\vec{p}) \quad \quad \quad
\textrm{subject to} \; \mathbf{A} \vec{p} = \vec{c} \; ,
\f}
where \f$\mathbf{A}\f$ is a \f$m\f$-by-\f$n\f$ matrix and \f$\vec{c}\f$ is a \f$m\f$-vector
with  \f$m \le n\f$.
In iterative methods the parameter vector \f$\vec{p}\f$ is expressed by
\f$\vec{p} = \vec{p}_k + \vec{d}\f$ with the correction \f$\vec{d}\f$
to \f$\vec{p}_k\f$ in the \f$k\f$-th iteration, satisfying the equation
\f{equation*}{
            \mathbf{A} \left( \vec{p}_k + \vec{d} \right) = \vec{c}
\f}
with \f$\vec{p}_{k+1} = \vec{p}_k +  \vec{d}\f$.

There are two methods for linear constrainst:
* in the *Lagrange multiplier method*  additional \f$m\f$ parameters
are introduced and the linear system of \f$n+m\f$ unknowns has
to be solved;
* by *elimination* the minimization problem with constraints
is transformed to
an unconstrained problem with \f$n-m\f$ unknowns. However the sparsity of
a matrix may be destroyed by the elimination.

The Lagrange method is used in **Millepede**.

\subsubsection sssec-lagrange The Lagrange multiplier method

In the Lagrange multiplier method  one additional parameter \f$\lambda\f$ is
introduced for each single constraint, resulting in an \f$m\f$-vector
\f$\vec{\lambda}\f$ of Lagrange multiplier. A term depending on \f$\vec{\lambda}\f$
and the constraints is added to the function  \f$F(\vec{p})\f$, resulting
in the Lagrange function
\f{equation*}
\mathcal{L}(\vec{p},\vec{\lambda}) = F(\vec{p}) + \vec{\lambda}
\left(  \mathbf{A} \vec{p} - \vec{c} \right)
\f}
Using as before a quadratic model for the function \f$F(\vec{p})\f$ and taking
derivatives w.r.t. the parameters \f$\vec{p}\f$ and the Lagrange multipliers
\f$\vec{\lambda}\f$, the two equations
\f{alignat*}{{2}
\mathbf{C} & \vec{d} +  \mathbf{A}^{\top} & \vec{\lambda} & = - \vec{g} \\
\mathbf{A} &  \vec{d}  &                                 & =
\vec{c} - \mathbf{A} \vec{p}_k
\f}
are obtained; the second of these equations is the constraint equation.
This system of two equations can be combined into one matrix equation

\anchor eq-lageq0 (10)
\f{equation*}{ \label{eq:lageq0}
\left(
    \begin{array}{ccc|c}
&           & &         \\
&    \mathbf{C} & & \mathbf{A}^{\top}  \\
&           & &         \\   \hline
&     \mathbf{A} & & \mathbf{0}
    \end{array}
    \right)
\left( \begin{array}{c}
~\\   \vec{d}  \\  ~\\  \hline \vec{\lambda}
        \end{array} \right)  =
\left( \begin{array}{c}
~\\  - \vec{g}  \\  ~\\  \hline   \vec{c} - \mathbf{A} \vec{p}_k
        \end{array} \right)  \; .
\f}
The matrix on the left hand side is still symmetric.
*Linear* least squares problems with *linear*
constraints can be solved directly, without iterations and without
the need for initial values of the parameters.

The matrix in equation \ref eq-lageq0 "(10)" is indefinite, with positive and
negative eigenvalues. A solution can be found even if the submatrix
\f$\mathbf{S}\f$ is singular, if the  matrix \f$\mathbf{A}\f$ of the constraints
supplies sufficient
information. Because of the different signs of the eigenvalues the
stationary solution is not a minimum of the function
\f$\mathcal{L}(\vec{p},\vec{\lambda})\f$.

\subsubsection sssec-feas Feasible parameters

<b> Feasible parameters. </b>
A particular value for the correction \f$\vec{d}\f$ can be  calculated by

\anchor eq-parsol (11)
\f{equation*}{ \label{eq:parsol}
\vec{d} =  \mathbf{A}^{\top} \left( \mathbf{A} \mathbf{A}^{\top} \right)^{-1}
\left(  \vec{c} - \mathbf{A} \vec{p}_k \right) \; ,
\f}
which is the *minimum-norm solution* of the constraint equation,
that is, the  solution of
\f{equation*}{
\min \; \left\| \mathbf{A} \left( \vec{p}_k + \Delta \vec{p} \right)
- \vec{c}  \right\|_2   \; ,
\f}
which is zero here.
The matrix \f$\mathbf{A}\f$ is a \f$m\f$-by-\f$n\f$ matrix for \f$m\f$ constraints and
the product \f$\mathbf{A} \mathbf{A}^{\top}\f$ is a square \f$m\f$-by-\f$m\f$ matrix,
which has to be inverted in equation \ref eq-parsol "(11)", which
allows to obtain a correction such that the linear constraints are satisfied.
Parameter vectors \f$\vec{p}\f$, satisfying the linear constraint
equations \f$\mathbf{A} \vec{p} = \vec{c}\f$, are called *feasible*.
If the vector \f$\vec{p}_k\f$ in the \f$k\f$-th iteration already
satisfies the linear constraint equations, then the correction
\f$\vec{d}\f$ has to have the property \f$\mathbf{A} \,\vec{d} = \vec{0}\f$.

\section sec_man The Manual

\subsection ssec_mp2 The programm package Millepede II

The second version of the program package with the name **Millepede**
(german: *Tausendf&uuml;ssler*) is based on the same mathematical
principle as the first version; global parameters are determined in
a simultaneous fit of global and local parameters. In the first version
the solution of the matrix equation for the global parameter corrections
was done by matrix inversion. This method is adequate
with respect to memory space and execution time for a number of global
parameters of the order of 1000. In actual applications e.g. for
the alignment of track detectors at the LHC storage ring at CERN
however the number of global parameters is much larger and may be
above \f$10^5\f$. Different solution methods are available in the
**Millepede** II version, which should allow the solution of problems
with a large number of global parameters with the
memory space, available in standard PCs, and with an execution time
of hours. The solution methods differ in the requirements of memory space
and execution time.

The structure of the second version **Millepede II** can be visualized
as the *decay* of the single program **Millepede** into two
parts, a part <b>%Mille</b> and a part **Pede** (Figure \ref fig-milped "1"):
\f{equation*}{
\rm{M}\small{ILLEPEDE} \; \Rightarrow \; \rm{M}\small{ILLE} \; + \; \rm{P}\small{EDE} \; .
\f}
The first part, <b>%Mille</b>, is a short subroutine, which is called in user
programs to write data files for **Millepede II**. The second part,
**Pede**, is a stand-alone program, which requires data files
and text files for the steering of the solution. The result is written
to text files.
\htmlonly <style>div.image img[src="fig_1.svg"]{width:600px;}</style> \endhtmlonly 
\image html fig_1.svg 


\subsubsection ssec_code Program code and makefile

The complete program is on a (\c tar)-file \c Mptwo.tgz, which
is expanded by the command

        tar -xzf Mptwo.tgz
into the actual directory. A \c Makefile for the
program *Pede* is included; it is invoked by the

        make
command.
There is a test mode **Pede** (see section \ref sssec-stalone), selected by

        ./pede -t
which can be used to test the program installation.

\anchor an-dynal The computation of the solution for a large number of parameter
requires large vector and matrix arrays. Most of the  memory space is required
in nearly all solution methods for a single matrix. The solution of
the corresponding matrix equation is done *in-space* for almost
all methods (no second matrix or matrix copy is needed).
The total space of all data arrays, used in **Pede**,  is defined
as a dimension parameter within the include file
\c dynal.inc by the statement

        PARAMETER       (MEGA=100 000 000) ! 100 mio words}
corresponding to the memory space allowed by
a 512 Mbyte memory. Larger memories allow an increase
of the dimension parameter in this statement in file \c dynal.inc.

Memory space above 2 Gbyte (?) can not be used with 32-bit systems, but
on 64-bit systems; a small change in the makefile to allow linking
with a big static array (see file \c dynal.inc) may be necessary
(see comment in makefile).

\subsubsection  sssec-mille1 Data collection in the user program with subroutine Mille

Data files are written within the user program by the subroutine
\c MILLE, which is available in Fortran and in C.
Data on the measurement and on derivatives with respect to local and
global parameters are written to a binary file.
The file or several files are
the input to the stand-alone program **Pede**, which performs
the fits and  determines the global parameters. The data required for
**Millepede** and the data collection in the user program
are discussed in detail in section \ref ssec-meapar.

\subsubsection sssec-stalone Solution with the stand-alone program Pede
The second part, **Pede**, is a *stand-alone program*,
which performs the fits and determines the global parameters.
It is written in Fortran.
Input to **Pede** are the binary (unformatted) data files, written using subroutine
<b>%Mille</b>, and text (formatted) files, which supply steering information and,
optionally, data about initial values and status of global parameters.
Different binary and text files can be combined.

\b Synopsis:

        pede [options] [main steering text file name]

The following options are implemented:
*  <tt>  -i   </tt>       interactive mode
*  <tt>  -t   </tt>       test mode
*  <tt>  -s   </tt>       subito option: stop after one data loop

<b> Option i.</b>
The interactive mode allows certain interactive steering, depending on
the selected method.

<b> Option t.</b>
In the test mode no user input files are required. This mode is recommended
to learn about the properties of **Millepede**.
Data files are generated by Monte Carlo simulation for a simple 200-parameter
problem, which is subsequently solved. The file generated are:
\c mp2str.txt (steering file), \c mp2con.txt (constraints) and
\c mp2tst.bin (datafile). The latter file is always (re-)created, the
two text files are created, if they do not exit. Thus one can edit these
files after a first test job in order to test different methods.

<b> Option s.</b>
In the *subito* mode, the **Pede**
program stops after the first data loop, irrespective of the options
selected in the steering text files.

<b> Text files.</b>
Text files should have either the characters \c xt or \c tx in
the 3-character filename-extension. At least one text file is necessary
(the default name  \c steer.txt is assumed, if no filename is given
in the command line), which specifies at least the
names of data files. The text files are described in detail in section
\ref ssec-textfiles.

**Pede** generates, depending on the selected method,
several output text files:

    millepede.log !  log-file for pede execution
    mpgparm.txt   !  global parameter results
    mpdebug.txt   !  debug output for selected events
    mpeigen.txt   !  selected eigenvectors

Existing files of the given names are renamed with a \f$\sim\f$ extension,
existing files with the \f$\sim\f$ extension  are removed.

\subsection ssec-meapar Measurements and parameters

\subsubsection sssec-mlp  Measurements and local parameters

Basic data elements are single measurements \f$y_i\f$; several single
measurements belong to a group of measurements, which can be called
a local-fit object. For example in
a track-based alignment based on tracks a track
is a local-fit object (see section \ref sssec-tracks for a detailed
discussion on the data for tracks in a track-based alignment).
A local-fit object is described by a linear or non-linear
function \f$f(x_i,\vec{q},\vec{p})\f$,
depending on a (small) number of local parameters \f$\vec{q}\f$ and in addition
on global parameters \f$\vec{p}\f$:

\anchor eq-measy (12)
\f{equation*}{ \label{eq:measy}
        y_i = f(x_i,\vec{q},\vec{p}) + \varepsilon_i  \; .
\f}
The parameters \f$\vec{q}\f$ are called the  local parameters,
valid for the specific group of measurements (local-fit object).
The quantity \f$\varepsilon\f$ is the measurement error, expected to have
(for ideal parameter values)
mean zero (i.e. the measurement is unbiased) and standard deviation
\f$\sigma_i\f$; often \f$\varepsilon\f$ will follow at least
The quantity \f$x_i\f$ is assumed to be the coordinate of the measured value
\f$y_i\f$, and is one of the arguments of the function \f$f(x_i,\vec{q},\vec{p})\f$.


The global parameters needed to compute the expectation
\f$f(x_i,\vec{q},\vec{p})\f$ are usually already quite accurate and only small
corrections habe to be determined. The convention is often to define
the vector \f$\vec{p}\f$ to be a correction with initial values zero. In this
case the final values of the global parameters will be small.

<b> The fit of a local-fit object.</b>
The local parameters
are usually determined in a least squares fit within the users code,
assuming fixed  global parameter values \f$\vec{p}\f$.
If the function
\f$f(x_i,\vec{q},\vec{p})\f$ depends *non-linearly* on the local parameters
\f$\vec{q}\f$, an iterative procedure is used, where the function
is linearized, i.e. the *first* derivatives of the function
\f$f(x_i,\vec{q},\vec{p})\f$ with respect to the local parameters \f$\vec{q}\f$
are calculated. Then the function is expressed as a linear function of local
parameter corrections \f$\vec{\Delta q}\f$ at some reference value \f$\vec{q}_k\f$:
\f{equation*}{
        f(x_i,\vec{q}_k +\vec{\Delta q},\vec{p})  =
        f(x_i,\vec{q}_k,\vec{p}) +
        \frac{\partial f}{\partial q_1} \Delta q_1
    + \frac{\partial f}{\partial q_2} \Delta q_2 + \ldots  \; ,
\f}
where the derivatives are calculated for \f$ \vec{q} \equiv \vec{q}_k\f$.
The corrections are determined by the linear least squares method,
a new reference value is obtained by
\f{equation*}{
        \vec{q}_{k+1} = \vec{q}_k + \vec{\Delta q}
\f}
and then, with \f$k\f$ increased by 1, this is repeated until convergence
is reached, i.e. until the corrections become essentially zero.
For each iteration a linear system of equations (normal equations of least
squares) has to be solved for the parameter corrections \f$\vec{\Delta q}\f$
with a right-hand side vector \f$\vec{b}\f$ with components
\f{equation*}{
        b_j= \sum_i \left( \frac{\partial f_i}{\partial q_j} \right)
\, \left(  y_i -  f(x_i,\vec{q}_k,\vec{p}_k)  \right) \frac{1}{\sigma_i^2}
\; ,
\f}
where the sum is over all measurements \f$y_i\f$ of the local-fit object. All
components \f$b_j\f$ become essential zero after convergence.

After convergence the equation \ref eq-measy "(12)" can be expressed with
the fitted parameters \f$\vec{q}\f$ in the form
\f{equation*}{
y_i = f(x_i,\vec{q},\vec{p}) +
        \frac{\partial f}{\partial q_1} \Delta q_1
    + \frac{\partial f}{\partial q_2} \Delta q_2 + \ldots
+ \varepsilon_i  \; .
\f}
The difference \f$z_i = y_i - f(x_i,\vec{q},\vec{p})\f$ is called the residual
and the equation can be expressed in terms of the
residual measurement \f$z_i\f$ as

\anchor eq-zdef (13)
\f{equation*}{ \label{eq:zdef}
z_i \equiv  y_i - f(x_i,\vec{q},\vec{p}) =
\left(  \frac{\partial f}{\partial q_1} \right) \Delta q_1
+  \left(  \frac{\partial f}{\partial q_2} \right)  \Delta q_2 + \ldots
+ \varepsilon_i  \; .
\f}
After convergence of the local fit all corrections
\f$\Delta q_1, \, \Delta q_2, \ldots\f$ are zero, and the residuals \f$z_i\f$ are small
and can be called the *measurement error*. As mentioned above the
model function \f$f(x_i,\vec{q},\vec{p})\f$ may be a linear or a non-linear
function of the local parameters. The fit procedure in the users code
may be a more advanced method than least squares. In track fits often
Kalman filter algorithms are applied to treat effects like multiple
scattering into account. But for any fit procedure the final result will
be *small* residuals \f$z_i\f$, and it should also be possible to
define the derivatives  \f${\partial f}/{\partial q_j}\f$ and
\f${\partial f}/{\partial p_{\ell}}\f$; these derivatives should express the
*change* of the residual \f$z_i\f$, if the local parameter \f$q_j\f$ or
the global parameter \f$p_{\ell}\f$ is changed by  \f$\Delta q_j\f$ or \f$\Delta p_{\ell}\f$.
In a track fit with strong multiple scattering the derivatives will become
rather small with increasing track length, because the information
content is reduced due to the multiple scattering.

Local fits are later (in **Pede**) repeated
in a simplified form; these fits need as information the difference \f$z_i\f$
and the derivatives, and correspond to the *last iteration* of the
local fit, because only small parameter changes ares involved.
Modified values of global parameters \f$\vec{p}\f$ during the global fit
result in a change of the value of \f$f(x_i,\vec{q},\vec{p})\f$ and
therefore of \f$z_i\f$, which can be calculated from the derivatives.
The derivatives with respect to the local parameters
allow to repeat the local fit and to calculate
corrections for the local parameters \f$\vec{q}\f$.

\subsubsection sssec-global Global parameters

Now the global parameters \f$\vec{p}\f$ are considered. The calculation
of the residuals \f$z_i = y_i - f(x_i,\vec{q},\vec{p})\f$ depends on the global
parameters \f$\vec{p}\f$ (valid for *all* local-fit objects).
This dependence can be considered in two ways. As expressed above,
the *parametrization* \f$f(x_i,\vec{q},\vec{p})\f$ depends on the global
parameters \f$\vec{p}\f$ and on corrections  \f$\vec{\Delta p}\f$ of the global
parameters (this case is assumed below in the sign of the derivatives with
respect to the global parameters). Technically equivalent is the case, where
the *measured* value \f$y_i\f$ is calculated from a raw measured value,
using global parameter values; in this case the residual could be
written in the form \f$z_i = y_i(\vec{p}) - f(x_i,\vec{q})\f$.
It is assumed that reasonable values are already assigned to the
global parameters \f$\vec{p}\f$ and the task is to find (small) corrections
\f$\vec{\Delta p}\f$ to these initial values.

Equation \ref eq-zdef "(13)" is extended to include corrections for
global parameters. Usually only few of the global parameters
influence a local-fit object.
A global parameter carries a label \f$\ell\f$,
which is used to identify the global parameter uniquely, and is
an arbitrary positive integer. Assuming that labels \f$\ell\f$ from a set
\f$\Omega\f$ contribute to a single measurement the extended equation
becomes

\anchor eq-zdefex (14)
\f{equation*}{ \label{eq:zdefex}
        z_i = y_i - f(x_i,\vec{q},\vec{p}) =
\sum_{j=1}^{\nu} \left( \frac{\partial f}{\partial q_j} \right)  \Delta q_j
+ \sum_{\ell \in \Omega}
\left( \frac{\partial f}{\partial p_{\ell}} \right) \Delta p_{\ell}  \; .
\f}
**Millepede** essentially performs a fit to all single measurements
for all groups of measurements simultaneously for all sets of local
parameters and for all
global parameters. The result of this simultaneous fit are optimal
corrections \f$\vec{\Delta p}\f$ for all global parameters; there are no
limitations in the number of local parameters sets. Within this global
fit the local fits (or better: the last iteration of the local fit) has
to be repeated for each local-fit object.

\anchor an-glolab <b> Global parameter labels.</b>
The label \f$\ell\f$, carried by a global parameter, is an arbitrary positive
(31-bit) integer (most-significant bit is zero), and
is used to identify the global parameter uniquely.
Arbitrary gaps between the numerical values of
labels values are allowed. It is recommended to design the label in a way
which allows to reconstruct the meaning of the global parameter from
the numerical value of the label \f$\ell\f$.

\subsubsection sssec-mille2 Writing data files with Mille

The user has to provide the following data for each single measurement:
\f{alignat*}{{2}
n_{lc} &= \;  \textrm{number of local parameters}
& \quad \quad
&\textrm{array}:  \; \left( \frac{\partial f}{\partial q_j}\right) \\
n_{gl} &= \;  \textrm{number of global parameters}
& \quad \quad
&\textrm{array}: \left( \frac{\partial f}{\partial p_{\ell}}\right) ;
    \; \textrm{label-array} \quad  \ell \\
z  &=  \;  \textrm{residual}  \quad
\left( \equiv y_i - f(x_i,\vec{q},\vec{p}) \right)
& \quad \quad
&\sigma =  \;  \textrm{standard deviation of the measurement}
\f}

These data are sufficient to compute the least squares normal equations
for the local and global fits.
Using calls of \c MILLE the measured data and the
derivatives with respect to the local and global parameters
are specified; they are collected in a buffer and
written as one record when one local-fit object is finished in the user
reconstruction code.

<b> Calls for Fortran version:</b> Data files should be written with the
Fortran \ref mille "MILLE" on the system used for **Pede**; otherwise the
internal file format could be different.
\verbatim
    CALL MILLE(NLC,DERLC,NGL,DERGL,LABEL,RMEAS,SIGMA)
\endverbatim

where
*       \c NLC \f$=\f$ number of local parameters \c DERLC
*       \c DERLC \f$=\f$ array \c DERLC(NLC) of derivatives
*       \c NGL \f$=\f$ number of global derivatives in this call
*       \c DERGL \f$=\f$ array \c DERLG(NGL) of derivatives
*       \c LABEL \f$=\f$ array \c LABEL(NGL) of labels
*       \c RMEAS \f$=\f$ measured value
*       \c SIGMA \f$=\f$ error of measured value (standard deviation)

After transmitting the data for all measured points the record
containing the data for one local-fit object is written.
following call
\verbatim
    CALL ENDLE
\endverbatim
<i> The buffer content is written to a file with the record.</i>

Alternatively the collected data for the local-fit object have to be
discarded,
if some reason for *not using* this local-fit object is found.
\verbatim
    CALL KILLE
\endverbatim
<i> The content of the buffer is reset, i.e. the data from preceeding
\c MILLE calls are removed.</i>

Additional floating-point and integer data (special data)
can be added to a local fit object by the call
\verbatim
    CALL MILLSP(NSPEC,FSPEC,ISPEC)
\endverbatim
where
* \c NSPEC \f$=\f$ number of special data \c DERLC
* \c FSPEC \f$=\f$ array \c FSPEC(NSPEC) of floating point data
* \c ISPEC \f$=\f$ array \c ISPEC(NSPEC) of integer data

The floating-point and integer arrays have the same length.
These special data are not yet used, but may be used in future options.

<b> Calls for C version:</b>
The C++-class Mille can be used to write C-binary files. The
constructor
\verbatim
Mille(const char *outFileName, bool asBinary = true, bool writeZero = false);
\endverbatim
takes at least one argument, defining the name of the output file. For
debugging purposes it is possible to give two further arguments: If
\c asBinary is false, a text file (not readable by pede) is written instead
of a binary output. If \c writeZero is true, derivatives that are zero
are not suppressed in the output as usual.

The member functions
\verbatim
void mille(int NLC, const float *derLc, int NGL, const float *derGl,
            const int *label, float rMeas, float sigma);
void special();
void end();
void kill();
\endverbatim
have to be called equivalently like the Fortran subroutines \c MILLE,
\c MILLSP, \c ENDLE} and \c KILLE.
To properly close the output file, the \c Mille object should be
deleted (or go out of scope) after the last processed record.

<b> Root version:</b> A \c root version is perhaps added in the future.

\subsubsection sssec-nonlin Non-linearities

Experience has shown that non-linearities in the function
\f$f(x_i,\vec{q},\vec{p})\f$ are usually small or at least not-dominating.
As mentioned before, the convention is often to define the start values
of the global parameters as zero and to expect only *small* final values.
Derivatives with respect to the local and global parameters can be
assumed to be constants within the range of corrections in the determination
of the global parameter corrections.
Then a single pass through the data, writing data files and calculating
global parameter corrections with  **Millepede** is sufficient.
If however corrections are large, they may require to re-calculate
the derivatives and another pass or several passes through the data
with re-writing data files
may be necessary. As a check of the whole procedure a second pass is
recommended.

\subsubsection sssec-tracks Application: Alignment of tracking detectors

In the alignment of tracking detectors the group of measurement,
the local-fit object,
may be the set of hits belonging to a single track. Two to five
parameters \f$\vec{q}\f$ may be needed to describe the dependence of the
measured values on the coordinate \f$x_i\f$.

In a real detector the trajectory of a charged particle does not follow
a simple parametrization because of effects like multiple scattering due
to the detector material. In practice often a Kalman filter method is
used without an explicit parametrization; instead the calculation proceeds
along the track taking into account e.g. multiple scattering effects.
Derivatives with respect to the (local) track parameters are not directly
available. Nevertheless derivatives as required in equation
\ref eq-zdefex "(14)" are still well-defined and can be calculated. Assuming
that the (local) track parameters refer to the starting point of the track,
for each point along the track the derivative
\f$(\partial f/\partial q_j)\f$ has to be equal to \f$\Delta z_i/\Delta q_j\f$, if the
local parameter \f$q_j\f$ if changed by \f$\Delta q_j\f$, and the
corresponding change of the
residual \f$z\f$ is \f$\Delta z_i\f$. This derivative could be calculated numerically.

For low-momentum tracks the multiple-scattering effects are large.
Thus the derivative above will tend to small values for measured points
with a lot of material between the measured point and the starting point of
the track. This shows that the contribution of low-momentum tracks
with large multiple-scattering effects to the precision of an alignment
is small.

\subsection ssec-textfiles Text files

In text files the steering information for the program execution
and further information on parameters and constraints is given.
The main text file is mandatory. Its default name is
\c steer.txt; if the name is different it has to be given
in the command
line of **Pede**. In the main text file the file names of data files
and optionally of further text files are given. Text files should have
an extension which contains the character pairs \c xt or \c tx.
In this section all options for the global fits and their keywords
are explained.
The computation of the global fit is described in detail
in the next section \ref ssec-globalfit; it may be necessary to
read this section \ref ssec-globalfit to get an understanding
of the various options.

\subsubsection sssec-dataform General data format of text files

Text files contain file names, numerical data, keywords (text) and
comment in free format, but there are certain rules.

<b> File names</b> for all text and steering file should be the leading
information in the main steering file, with one file name per line,
with correct characters, since the file names are used to open the files.

<b> Numerical </b> can be given with or without decimal point,
optionally with exponent field. The numerical data

    13234       13234.0      13.234E+3

are all identical.

<b>Comments</b> can be given in every line, preceeded by the ! sign;
the text after the ! is ignored. Lines with * or ! in the first
columm are considered as comment lines. Blank lines are ignored.

<b>Keywords</b> are necessary for certain data; they can be given in upper or
lower case characters. Keywords with few typo errors may also be
recognized correctly. The program will stop if a keyword is not
recognized.

Below is an example for a steering file:
<small>\verbatim
Fortranfiles
!/home/albert/filealign/lhcrun1.11    ! data from first test run
/home/albert/filealign/lhcrun2.11    ! data from second run
/home/albert/filealign/cosmics.bin   ! cosmics
/home/albert/detalign/mydetector.txt ! steering file from previous log file
/home/albert/detalign/myconstr.txt   ! test constraints

constraint 0.14       ! numerical value of r
713 1.0               ! pair of parameter label and numerical factor
719 0.5  720 0.5      ! two pairs

CONSTRAINTS 1.2
112 1.0
113 -1.0
114 0.5
116 -0.5

Parameter
201 0.0   -1.0
202 1.0   -1.0
204 1.23 0.020

method inversion 5 0.1
end
\endverbatim  </small>

\subsubsection sssec-fileinf File information

Names of data and text files are given in single text lines.
The file-name extension is used to distinguish text files (extension
containing \c tx or \c xt) and binary data files.
Data files may be Fortran or C files; by default C files
are assumed. Fortran files have to be preceded by the textline

    Fortranfiles

and C files may be preceded by the textline

    Cfiles

Mixing of C- and Fortran-files is allowed.

\subsubsection sssec-parinf Parameter information

By default all global parameters appearing in data files with their
labels are assumed to be
variable parameters with initial value zero, and used in the fit.
Initial values different from zero and the so-called presigma
(see below) may be
assigned to global parameters, identified by the parameter label.
The optional parameter information has to be provided in the form below,
starting with a line with the keyword \c Parameter:

    Parameter
    label   initial_value   presigma
        ...
    label   initial_value   presigma

The three numbers \a label, \a initial_value and \a presigma
in a textline may be followed by further numbers, which are
ignored.
\note Note that the result file has the same format and
may be used as input file, if a program execution should be continued.

All parameters not given in this file, but present in the data files
are assumed to be variable, with initial value of zero.

<b>Pre-sigma.</b>
The pre-sigma \f$s_{\ell}\f$ defines the status of the global parameter:

\f$\boldsymbol{s_{\ell} > 0}\f$:
The parameter is *variable* with the given initial value.
A term \f$1/s_{\ell}^2\f$ is added to the diagonal matrix element of the
global parameter to stabilize a perhaps poorly defined parameter.
This addition should *not* bias the fitted parameter value, if a
sufficient number of iterations is performed; it may however bias the
calculated error of the parameter (this is available only for the
matrix inversion method).

\f$\boldsymbol{s_{\ell} = 0}\f$:
The pre-sigma is zero; the parameter is *variable* with the
given initial value.

\f$\boldsymbol{s_{\ell} < 0}\f$:
The pre-sigma is negative. The parameter is defined as \b fixed; the
initial value of the parameter is used in the fits.

The lines below show examples, with  the numerical value of the label,
the initial value, and the pre-sigma:
<small>
\verbatim
11 0.01232  0.0    ! parameter variable, non-zero initial value
12 0.0      0.0    ! parameter variable, initial value zero
20 0.00232  0.0300 ! parameter variable with initial value 0.00232
30 0.00111  -1.0   ! parameter fixed, but non-zero parameter value
\endverbatim
</small>

<b> Result file.</b>
At the end of the **Pede** program a result file \c millepede.res
with the result
is written. In the first three columns it carries the label, the fitted
parameter value (initial value + correction) and the pre-sigma (default is zero). This file can,
perhaps after editing, be used for a repeated **Pede** job. Column 4 contains the correction and
in the case of solution by inversion or diagonalization column 5 its error
and optionally column 6 the global correlation.

\subsubsection sssec_consinf Constraint information

*Equality constraints*
allow to fix certain undefined or weakly defined
linear combinations of global parameters by the addition of
equations of the form
\f{equation*}{
c = \sum_{\ell \in \Omega} f_{\ell} \cdot p_{\ell}
\f}
The constraint value \f$c\f$ is usually zero. The format is:

    Constraint   value
    label        factor
        ...
    label        factor

where \a value, \a label and \a factor
are numerical values. Note that no error of the value is given.
The equality constraints are, within the numerical accuracy, exactly
observed in a global fit. Each constraint adds another Lagrange multiplier
to the problem, and needs the corresponding memory space in the matrix.

Instead of the keyword \c Constraint the keyword \c Wconstraint
can be given (<tt>this option is not yet implemented!</tt>).
In this case the factor for each global parameter given
in the text file is, in addition,
multiplied by the weight \f$W_{\ell}\f$ of the corresponding global parameter
in the complete fit:
\f{equation*}{
c = \sum_{\ell \in \Omega} f_{\ell} \cdot W_{\ell} \cdot  p_{\ell}
\f}
Thus global parameters with a large weight in the fit
get also a large weight in the constraint.

Mathematically the effect of a constraint does not change, if the
constraint equation is multiplied by an arbitrary number. Because of
round-off errors however it is recommended to have the constraint equation
on a equal accuracy level, perhaps by a scale factor for
constraint equations.

\subsubsection sssec_gpm Global parameter measurements

*Measurements* with measurement *value* and measurement *error* for linear
combinations of global parameters in the form
\f{equation*}{
y = \sum_{\ell \in \Omega} f_{\ell} \cdot p_{\ell} + \epsilon
\f}
can be given, where \f$\epsilon\f$ is the measurement error, assumed to
follow a distribution \f$N(0,\sigma^2)\f$, i.e. zero bias and standard deviation
\f$\sigma\f$.
The format is:

    Measurement    value    sigma
    label    factor
        ...
    label    factor

where \a value, \a sigma, \a label and \a factor
are numerical values. Each measurement
\a value \f$\pm\f$ \a sigma (standard deviation) given in the text files
contributes to the overall objective function in the global fit.
At present no outlier tests are made.

The global parameter *measurement*
and the *constraint* information from the
previous section differ, in the definition, only in the
additional information on the standard deviation
given for the measurement. They differ however in the mathematical
method: constraints add another parameter (Lagrange multiplier);
measurements contribute to the parameter part of the matrix, and may
require no extra space; however a sparse matrix may become less
*sparse*, if matrix elements are added, which are empty otherwise.

\subsubsection  sssec-methodsel Solution method selection

<b> Methods.</b>
One out of several solution methods can be selected. The methods have
different execution time and memory space requirement, and may also have
a slighly different accuracy. The statement to select a method are:

    method inversion         number1 number2
    method diagonalization   number1 number2
    method fullGMRES         number1 number2
    method sparseGMRES       number1 number2
    method cholesky          number1 number2
    method bandcholesky      number1 number2
    method HIP               number1 number2

The two numbers are:
* \a number1 = number of iterations
* \a number2 = limit for \f$\Delta F\f$ (convergence recognition).

For preconditioning in the GMRES methods a band matrix is used.
The width of the variable-band matrix is defined by

    bandwidth     number

where \a number is the numerical value of the semi-bandwidth of
a band matrix. The method called \c bandcholesky
uses only the band matrix
(reduced memory space requirement), and needs more iterations.
The different methods are described below.

<b>Iterations and convergence.</b>
The mathematical method in principle allows to solve the problem in a
single step. However due to potential inaccuracies in the solution of
the large linear system and due to a required outlier treatment certain
internal iterations may be necessary.
Thus
the methods work in iterations, and each iteration requires one or several
loops with
evaluation of the objective function and the first-derivative vector
of the objective function (gradient), which requires reading all data
and local
fits of the data. The first loop is called iteration 0, and (only) in this
loop in addition the second-derivative matrix is evaluated, which takes
more cpu time. This first loop usually brings a large reduction of the
value of the objective function, and all following loops will bring
small additional improvements.
In all following loops the same second-derivative matrix
is used (in order to reduce cpu time), which should be a good approximation.

A single loop may be sufficient, if there are no outlier problems.
Using the command line option \f$\bf -s\f$ (*subito*) the program
executes only a single loop, irrespective of the options required
in the text files.
The treatment of outliers turns a linear problem into a non-linear problem
and this requires iterations. Also precision problems, with
rounding errors in the case of a large number of global parameters,
may require iterations. Each iteration \f$1, \, 2, \ldots\f$ is a line search
along the calculated search direction, using the *strong Wolfe
conditions*, which force a *sufficient* decrease of the function
value and the gradient. Often an iteration requires only one loop.
If, in the outlier treatment, a cut is changed, or if the
precision of the constraints is insufficient,  one additional
loop may be necessary.

At present the iterations end, when the number specified with the method
is reached. For each loop, the expected objective-function decrease and the
actual decrease are compared. If the two values are below the limit
for \f$\Delta F\f$ specified with the method, the program will end earlier.

<b>Memory space requirements.</b>
The most important consumer of space is the symmetric matrix of the
normal equations. This matrix has a parameter part, which requires
for example
for full storage \f$ n(n+1)/2\f$ words, and for sparse storage
\f$ n + q n(n-1)/2\f$ words, for \f$n=\f$ number of global parameters and
\f$q =\f$ fraction of non-zero off-diagonal
elements. In addition the matrix has a constraint part, corresponding
to the Lagrange multipliers. Formulae to calculate the number of
(double precision) elements of the matrix are given in table
\ref tab-space "1".

The GMRES method can be used with a full or with a sparse matrix.
Without preconditioning the GMRES method may be slow or may even fail
to get an acceptable solution. Preconditioning is recommended and
usually speeds up the GMRES method considerably. The band matrix required
for preconditioning also takes memory space. The parameter part requires
(only) \f$n \cdot m\f$ words, where \f$m =\f$ semibandwidth specified with the
method. A small value like \f$m =6\f$ is usually sufficient. The constraint
part of the band matrix however requires up to
\f$n \cdot n_C + n_C (n_C+1)\f$ words, and this can be a large number for
large number \f$n_C\f$ of constraint equations.

Another matrix with \f$n_C (n_C +1)/2\f$ words is required for the constraints,
and this matrix is used to improve the precision of the constraints.
Several vectors of length \f$n\f$ and \f$n_C\f$ are used in addition, but this
space is only proportional to \f$n\f$ and \f$n_c\f$.

The diagonalization method requires more space:
in addition to the symmetric matrix
another matrix with \f$(n+n_C)^2\f$ words is used for the eigenvectors;
since diagonalization
is also slower, the method can only be used for smaller problems.

\anchor tab-space Table 1
\image html table-space.png " "
\latexonly
\begin{table}
\begin{center}\begin{tabular}{ll} \toprule
method                        & memory space requirement \\ \midrule
inversion, fullGMRES, cholesky & $ (n^2+n)/2 + n n_C +
                                (n_C^2 +n_C)/2$ \\
diagonalization                & $         n + n(n-1)/2 + n^2$ \\
sparseGMRES                    & $ n + q n(n-1)/2
                                n_{cf}   $ \\
bandcholesky, preconditioning     & $ nm
                            + n n_C + (n_C^2 +n_C)/2  $ \\
HIP (no constraints)           & $nm$ \\  \bottomrule
\end{tabular} \end{center}
\caption*{ \label{tab:space}
Table 1: Matrix space requirements for the different methods, with:
$q =$ fraction of non-zero off-diagonal elements, $m=$ bandwidth of
variable band matrix, $n_C =$ number of constraints, $n_{cf}=$ number of
constraint-factors.}
\end{table}
\endlatexonly

<b>Comparison of the methods.</b> The methods have a different
computing-time dependence on the number of parameters \f$n\f$. The methods
\c diagonalization, \c inversion and \c cholesky have a
computing time proportional to the third power of n, and can therefore
be used only up to \f$n=\f$ 1000 or perhaps up to \f$n=\f$ 5000. Especially the
diagonalization is slow, but provides a detailed information on
weakly- and well-defined linear combinations of global parameters.

The GMRES methods has a weaker dependence on \f$n\f$, and allows much larger
\f$n\f$-values; larger values of \f$n\f$ of course require a lot of space and
the sparse matrix mode should be used. The bandcholesky
method can be used for preconditioning the GMRES method (recommended),
but also as a stand-alone method; its cpu time depends only linearly
on \f$n\f$, but the matrix is only an approximation and usually many
iterations are required. For all methods the constraints increase the
time and space consumption; this remains modest if the number of
constraints \f$n_C\f$ is modest, for example \f$n_C \le 100\f$, but may be
large for \f$n_C =\f$ several thousand.

\subsection sssec-methods Solution methods

\anchor an-inv <b>Inversion.</b> The computing time for inversion is roughly
\f$2 \times 10^{-8} \cdot n^3\f$ seconds (on a standard PC).
For \f$n \approx 1000\f$ the inversion time
of \f$ \approx 20\f$ seconds
is still acceptable, but for \f$n \approx 10000\f$ the inversion time would be
already more than 5 hours.
The advantage of the inversion is the availability of *all*
parameter errors and global correlations.

<b>Cholesky decomposition.</b> The Cholesky decomposition (under test)
is an alternative
to inversion; eventually this method is faster and/or numerically more
stable than inversion. At present parameter errors and global correlations
are not available.

\anchor an-diag <b>Diagonalization.</b> The computing time for diagonalization is
roughly a factor 10 larger than for inversion. Space for the square
(non-symmetric) transformation matrix of eigenvectors is necessary.
Parameter errors and
global correlations are calculated for all parameters. The main advantage
of the diagonalization method is the availability of the eigenvalues.
Small positive eigenvalues of the matrix
correspond to linear combinations of parameters
which are only weakly defined and it will become possible
to interactively
remove the weakly defined linear combinations. This removal could be an
alternative to certain constraints. Constraints are also possible and they
should correspond to *negative* eigenvalues.

\anchor an-gmres <b>Generalized minimization of residuals (GMRES).</b> The GMRES method
is a fast iterative solution method for full and sparse matrices.
Especially for large dimensions with \f$n \gg\f$ 1000 it should be much faster
than inversion, but of similar accuracy.
Although the solution is fast, the building of the sparse matrix
may require a larger cpu time.
At present there is a limit
for the number of internal GMRES-iterations of 2000. For a matrix with a
bad condition number this number of 2000 internal GMRES-iterations will
often be reached and this is often also an indication of a bad and probably
inaccurate solution. An improvement is possible by preconditioning,
selected if GMRES is selected together with a bandwidth parameter
for a band matrix; an approximate solution is determined by solving the
matrix equation with a band matrix within the
GMRES method, which improves the eigenvalue spectrum of the matrix and
often speeds up the method considerably. Experience shows that
a small bandwidth of e.g. 6 is sufficient.
In interactive mode parameter errors for selected parameters can be
determined.

<b>Variable band matrix.</b> Systems of equations with a band matrix
can be solved by the Cholesky decomposition. Often the matrix element
around the diagonal are the essential matrix elements and in these
cases the matrix can be approximated by a band matrix of small
width, with a small memory requirement. The computing time is also
small; it is linear with the number of parameters.
However because
the band matrix is only the approximate matrix often many iterations
are necessary to get a good solution. The definition of the band width
is necessary.
The band width
is variable and thus the constraints equations can be treated
completely. This means that the constraints are observed even in this
approximate solution.

\subsubsection sssec-outlierdeb  Outlier treatment and debugging

Cases with very large residuals within the data can distort the result of
the least squares fit. Reason for outliers may be selection mistakes or
statistical fluctuations. A good outlier treatment may be neccessary to get
accurate results.
The problem of outlier treatment is complicated
by the fact that initially the global parameters may be far from optimal
and therefore large deviation may occur even for correct data before the
global parameter determination.
There are options for the complete removal of bad cases and for the
down-weighting of bad data. Cases with a
        <b>huge \f$\chi^2\f$ are automatically removed</b>
in every iteration. The options to remove large \f$\chi^2\f$ cases
and down-weighting are *not* done  in the first iteration. The
options are:

    chisqcut              number1 number2
    outlierdownweighting  number
    dwfractioncut         number
    printrecord           number1 number2

\anchor an-chisq <b>Chisquare cut.</b> With the keyword \c chisqcut two numbers can
be given. Basis of \f$\chi^2\f$ rejection is the \f$\chi^2\f$-value
corresponding to 3 standard deviations (and to a probability of
\f$0.27 \%\f$).
or one degree of freedom this \f$\chi^2\f$-value
is 9.0, and for 10 degrees of freedom  the \f$\chi^2\f$-value is 26.9. The first
number given with the keyword \c chisqcut is a factor for the
\f$\chi^2\f$-cut value, to be used in the first iteration, and the second
number is the factor for the second iteration. In subsequent iterations
the factor is reduced by the square root, with values below 1.5 replaced
by 1. For example the statement

    chisqcut 5.0 2.5

means: the cut factor is \f$5 \times \chi^2\f$-value corresponding to three
standard deviations, \f$2.5 \times \chi^2\f$-value for the second iteration,
\f$1.58 \times \chi^2\f$-value for the third iteration and
\f$1 \times \chi^2\f$-value for subsequent iterations. Thus in the first
iteration the cut is
\f$5 \times 26.9 = 134.5\f$ for 10 degrees of freedom, and \f$26.9\f$ after
the third iteration.

\anchor an-downw <b>Outlier downweighting.</b>
Outlier down-weighting for single data
values requires repeated local fits with \f$>1\f$ iterations in the local fit.
The number given with the keyword \c outlierdownweighting is the number
of iterations.
In down-weighting the weight of the data
point, which by default is \f$1/\sigma^2\f$, is reduced by an extra factor
\f$\omega_i < 1\f$ according to the residual (see below).
For \f$n\f$ data points the total sum \f$S_f = \sum_i \omega_i\f$ of the extra
factors is \f$\le n\f$.
The ratio \f$(n - S_f)/n\f$ is called the down-weight fraction, and a
histogram of this ratio is accumulated in the second function evaluation.
The first iteration of the local fit is done *without*
down-weighting, iterations 2 and 3 use the M-estimation method with
the Huber function. Subsequent iterations, if requested, use the
M-estimation method with the
Cauchy function for down-weighting, where the influence of very large
deviations is further reduced. For example the statement

    outlierdownweighting 4

means: iteration 1 of the local fit without down-weighting, iterations
2 and 3 with Huber function down-weighting and iteration 4 with
Cauchy function down-weighting.

\anchor an-dwcut <b>Downweighting fraction cut.</b> Cases with a very large
fraction of down-weighted measurements are very likely wrong data and
should be rejected. The cut value should be determined
from the histogram showing the down-weight fraction. Typical values
are 0.1 for weights from the Huber function
(\c outlierdownweighting
\f$\le 3\f$) and 0.2 for weights from the Cauchy function
(\c outlierdownweighting \f$> 3\f$). For example the statement

    dwfractioncut 0.2

means to reject all cases with a down-weight fraction \f$\ge 0.2\f$.

\anchor an-recpri <b>Record printout.</b> For debugging many details of single
records and the local fits are printed. Two record number  <i>number1  number2</i>
can be given with the
keyword \c printrecord; printout is done in iteration 1 and 3.
For numbers given as negative,
the records with  extreme deviations are selected
in iteration 2, and printed in iteration 3. If the first number is given
as \f$-1\f$, the record with the largest single residual is selected. If the
second number is given as \f$-1\f$, the record with the largest value of
\f$\chi^2/N_{df}\f$ is selected. For example the statement

    printrecord 4321 -1

means: the record 4321 is printed in iterations 1 and 3, and the
record with the largest value of \f$\chi^2/N_{df}\f$ is selected in iteration
2 and printed in iteration 3.

\subsubsection sssec_further Further options

There are miscellaneous further options which can be selected within the text
files:

    subito
    entries           number
    nofeasiblestart
    wolfe             C1 C2
    histprint
    end

The option \c subito is also selected by \b -s in the command line.
The program will end regularly, with normal output of the result, already
after the first function evaluation (iteration 0).

\anchor an-entries The number, given with the keyword \c entries,
is the minimum number of
data for a global parameter. Global parameters which have a
number of entries \f$=\f$ number of measured points
connected to it in the local fit-objects, smaller then the given number
are ignored. This option is used to suppress global parameters
with only few data, which are therefore rather inaccurate, and would
spoil the condition of the linear system.

\anchor an-nofeas By default the initial global parameter values are corrected
to follow all the constraints. A check for the validity of the constraints is
repeated in every iteration, and eventually another correction is made
at the end of an iteration.
With the keyword \c nofeasiblestart the
parameters are *not* made feasible (respecting the constraints)
at the start.

\anchor an-wolfe The constants \a C1 and \a C2 are the line search
constants for the strong Wolfe conditions; default values are the
standard values \f$C_1 = 10^{-4}\f$ and \f$C_2 = 0.9\f$.

\anchor an-histpr Histograms accumulated during the program are written to the
textfile \c millepede.his and can be read after the job execution.
If selected by the keyword \c histprint the histograms
accumulated during the program are also printed.

The reading of a textfile ends, when the keyword \c end is
encountered. Textlines after the line with \c end are not read.

\subsection ssec-globalfit Computation of the global fit

\subsubsection sssec-memman Memory managament

The solution of the optimization problem for a large number of
parameters requires one large memory with many arrays,
required to store vectors and matrices of large size corresponding to the
number of parameters and constraints. **Millepede II** uses a large
array in a common. This large array is dynamically divided into
so-called subarrays, which are created, enlarged and moved,
and removed according to the requirements.
The space limitation is thus given for almost all
subproblems only by the total space requirement. The largest space
is required for the symmetric matrix of the normal least squares
equations. This matrix can be a full matrix or for certain
problems a sparse matrix. As an approximation a band matrix with a small
band width can be used, if there is not enough space for the full or
sparse matrix.

\subsubsection sssec-init Initialization

After initialization of the memory management the command line options are
read and analysed. Then the main text file is analysed, the names
of other text files and the data files are recognized and stored in a table.
Default value are assigned to the various variables; then all text files
are read and the requested options are interpreted. The data for
the definition for parameters, constraints and measurement are read and
stored in subarrays.

\subsubsection sssec-loop1 First data loop

In subroutine \c LOOP1 tables for the variable and
fixed global parameters and translation tables are defined. The table
of global parameter labels is first filled with the global parameters,
appearing in the text files. All data files are read and all global
parameter labels from the records are included in the list.

There are three integers to characterize a global parameter.

<b>Global parameter label = \c ITGBL:</b> this is a positive integer,
from the full range of 32-bit integers.
\f$1 \ldots 2147483647= 2^{31} -1\f$.

<b>Global parameter index = \c ITGBL:</b>
A translation
table is constructed to translate the global parameter label
\c ITGBL to a global parameter index \c ITGBI, with a range
\f$1 \ldots\f$ \c NTGB, where \c NTGB is the total number of global
parameters  (variable and fixed parameters).

<b>Variable-parameter index = \c IVGBI:</b> the global parameters which
are *variable* carry a variable-parameter index, with a range
\f$1 \ldots\f$ \c NVGB, where \c NVGB is the number of variable global
parameters. Parameter vectors in the mathematical computation are
vectors containing only the variable parameters.

Function calls are used to translate the global parameter label \c ITGBL:
\f{alignat*}{{2}
\textrm{global parameter index}
&\leftarrow
\textrm{global parameter label}
        & \quad \quad \quad
\texttt{ITGBI} &= \texttt{INONE(ITGBL)}  \\
\textrm{variable parameter index}
&\leftarrow
\textrm{global parameter label}
        & \quad \quad \quad
\texttt{IVGBI} &= \texttt{INSEC(ITGBL)}  \; .
\f}
The translation is done with a hash-index table.
The translation from one parameter index
to another one and back to the parameter label
is done by the following statement functions:
\f{alignat*}{{2}
\textrm{global parameter label}
&\leftarrow
\textrm{global parameter index}
        & \quad \quad \quad
    \texttt{ITGBL} &= \texttt{JTGBL(ITGBI)}   \\
\textrm{variable-parameter index}
&\leftarrow
\textrm{global parameter index}
    & \quad \quad \quad
    \texttt{IVGBI} &= \texttt{JVGBI(ITGBI)}  \\
\textrm{global parameter index}
&\leftarrow
\textrm{variable-parameter index}
        & \quad \quad \quad
    \texttt{ITGBI} &= \texttt{JTGBI(IVGBI)}  \; ,
\f}
which use simple fixed-length tables.

\subsubsection sssec-loop2 Second data loop

In subroutine  \c LOOP2
the subarray dimensions of the various vectors and matrices are determined.
All data files are read and for each local-fit object the number
of local and of global parameters is determined; the maximum values
of the numbers are used to define the subarray dimensions for the
arrays to be used in local fits, and for the contributions to the global
fit.

The sparse matrix storage requires to establish the pointer structure
of the sparse matrix. During reading the data files the list of global
parameters used in each local-fit object is determined; an entry is made
in a table for each pair of global parameters, which later corresponds
to an off-diagonal matrix element. A hash-index table is used for the
search operation.
For sparse and full matrix storage the requirements of the constraint
equations and Lagrange multipliers has to be included.
For the sparse matrix a pointer structure is
build, which allows a fast *matrix*\f$\times\f$*vector* product.
At the end of the subroutine all subarrays required for the determination
of the solution are prepared.

\subsubsection sssec-solmetover Solution method overview

The section gives an short overview over the solution method; more detailed
explanation is given in the subsequent sections.

In principle the solution can be determined in a single step; this
mode can be selected by the keyword \c subito or option \b -s.
There are however two reasons to perform iterations, which improve the
result:
* Due to the large size or the method  the one-step solution may be
affected by rounding errors and may not be precise. Experience shows that
the overall value objective-function value can be reduce using
more than one step, although the decrease is sometimes rather small;
* Outliers may be important; they add a *non-linear component to the
otherwise linear problem and this requires iterations and repeated
evaluation of the objective function; each function evaluation
requires to read again all data files and to
repeat the local fits. In a comparison of Millepede I and II
results initially certain differences were observed; only after a
careful outlier treatment these differences became small or disappeared.

The solution  determined in iterations is explained.
The starting iteration, with iteration number 0, is the most important and
time-consuming one. The data files are read, and for each case a local fit
is made. The matrix of the normal equations is formed only in this starting
iteration. Depending on the selected method the matrix is inverted,
diagonalized or decomposed, and the resulting matrix and decomposition
is used in all later data-file loops, which take less time compared to the
first data-file loop.
Thus the data-file loop in iteration number 0 may be rather time
consuming, especially for the case of a sparse matrix, where the index
determination takes some time. The right-hand-side vector is formed in each
data-file loop.
A correction step in global parameter space is calculated
at the end of the data-file loop in iteration
number 0. Often the
precision is already sufficient and no further data-file loops are
necessary. In the subito mode (keyword \c subito or option \b -s)
the correction is added to the global parameter values and the program stops.

Iterations with more data-file loops are recommended, to check the result,
eventually to improve the result and to treat outliers. In each iterations
a so-called line search is done, where the overall value of the objective
functions is optimized along the correction-step direction
*sufficiently*, using the strong Wolfe criterion. Often a single
step is already sufficient for an iteration. The sample of local-fit
objects may change during the iterations becuse of changing cut values;
this may require extra function and derivative calculations.

\subsubsection sssec-loopn Data file loops during iterations

In subroutine \c LOOPN all data files are read. For each local-fit
object the local fit is performed, based on the actual global parameter
corrections, and eventually including downweighting of outliers.
Using the results of the local fit the value of the objective function
\f$F\f$, the vector of first derivatives (gradient) and, in the first
data file loop, the matrix of second derivatives is calculated. For a linear
fit the matrix of second derivatives will not change during the iterations;
the outlier treatment will change the matrix, but usually the changes are
small und the matrix collected in the first data loop is at least a good
and sufficiently accurate approximation. Thus data loops after the first
one are faster; depending on the method the solution of the matrix equation
will be faster after the first data loop.

\paragraph par-locfitv The local fit

The number of local parameters of a local-fit object is usually small;
typical values are between 2 and 5. Since the problem is linear, a single
step is sufficient in a local fit with minimization of the sum
\f{equation*}{
\tfrac{1}{2} \sum_i \left( \frac{y_i - f(x_i,\vec{q},\vec{p})}{\sigma_i}
\right)^2  \; ,
\f}
unless outlier downweighting is required, which may require a few
iterations.
In a loop
over the local measurements the matrix and the right-hand side of the
least squares normal equations are summed. For each single measured value,
the residual measurement \f$z_i\f$
\f{equation*}{
z_i \equiv  y_i - f(x_i,\vec{q},\vec{p})
\f}
(see equation \ref eq-zdef "(13)")
has to be corrected for the actual global parameters corrections
\f$\vec{\Delta p}\f$ using the first global parameter derivatives (see equation
\ref eq-zdefex "(14)").
The corrected residual \f$z_i'\f$ is then used in the accumulation of the
matrix and vector:
\f{equation*}{
\Gamma_{jk} = \sum_i \left( \frac{\partial f_i}{\partial q_j} \right)
        \left( \frac{\partial f_i}{\partial q_k} \right)
            \frac{1}{\sigma_i^2}
            \quad \quad \quad \quad
        \beta_j=  \sum_i \left( \frac{\partial f_i}{\partial q_j} \right)
\,  \frac{z_i'}{\sigma_i^2}  \; .
\f}
Corrections \f$\vec{\Delta q}\f$ are determined by the solution of the
matrix equation
\f{equation*}{
            \vec{\Gamma} \vec{\Delta q} = - \vec{\beta} \; ,
\f}
which is determined using matrix inversion
\f$\vec{\Delta q} =  - \vec{\Gamma}^{-1} \vec{b}\f$, because the inverse matrix
\f$\vec{\Gamma}^{-1}\f$ (the covariance matrix)
is necessary for the contribution to
the matrix of the global normal equations. The residual \f$z_i'\f$ are then
corrected for the local parameter corrections:
\f{equation*}{
z_i^{''} = z_i^{'} - \sum_j \frac{\partial f_i}{\partial q_j} \Delta q_j
\f}
and the new residuals \f$z_i^{''}\f$ are used to calculate the
\f$\chi^2\f$ value \f$S\f$ of the local fit
\f{equation*}{
S = \sum_i \left( \frac{z_i^{''}}{\sigma_i} \right)^2  \; ,
\f}
which should follow a \f$\chi^2\f$. The number of local parameters of a local-fit
object is usually small;
typical values are between 2 and 5. Since the problem is linear, a single
step is sufficient in a local fit with minimization of the sum
\f{equation*}{
\tfrac{1}{2} \sum_i \left( \frac{y_i - f(x_i,\vec{q},\vec{p})}{\sigma_i}
\right)^2  \; ,
\f}
unless outlier downweighting is required, which may require a few
iterations. In a loop
over the local measurements the matrix and the right-hand side of the
least squares normal equations are summed. For each single measured
distribution with the given number of
degrees of freedom \f$n_{\textrm{df}} =\f$ number of measurements minus
number of parameters. A fraction of 0.27 \% or one out of 370
of the *correct* cases should have a deviation which corresponds
to 3 standard devations in the \f$n_{\textrm{df}} =1\f$ case.

\anchor localfit-rejection
Very badly fitting cases should be rejected before using them for
the global parameter fit.  Basis of the rejection is the comparison
of the sum \f$S\f$ and the *0.27 \% value*
\f$\chi^2_{\textrm{cut}}\f$   of the  \f$\chi^2_{n_{\textrm{df}}}\f$
distribution. If the value \f$S\f$ exceeds \f$\chi^2_{\textrm{cut}}\f$ by more than
a factor of 50, then the sum is called *huge* and the local-fit
object is rejected. Cases with \f$n_{\textrm{df}} =0\f$ are rejected too,
because they do no allow any check. These rejections are always made,
even if the user has not requested any outlier rejection.

With keyword \c chisqcut the user can require a further rejection.
The first number given with the keyword is used during iteration 0;
the local-fit object is rejected,
if the value \f$S\f$ is larger than this number times
\f$\chi^2_{\textrm{cut}}\f$.
Often this number has to rather high, e.g. 12.0;
the initial value of \f$S\f$ may be rather high before the first
determination of global parameters, even for correct data, and, if possible,
no correct data should be rejected by the cut. The second number given with
the keyword is used during iteration 1 and can be  smaller than the first
number, because the first improvement of the local parameters is large.
In further
iterations the number is reduced by the sqrt-function. Values below
1.5 are replaced by 1, and thus the rejection is finally be done
with a 0.27 \%-rejection probability for correct data.

All not rejected local-fit objects contribute to the global fit.
The function value \f$F\f$
is the sum of the \f$S\f$-values for the local-fit objects:
\f{equation*}{
            F = \sum_{l} S_{l}
\quad \quad \quad \textrm{sum with index $l$ over all local fit-objects}
\; .
\f}
Here even \f$S\f$-values from local-fit objects, which are rejected
by cuts are included, with \f$S_{l} = \f$ cut value, in order
to keep the sum \f$F\f$ meaningful (otherwise the rejection of many local-fit
object could simulate an improvement -- which is not the case).

\paragraph par-conglo Contributions to the global parameter fit

After an accepted local fit the contributions to the normal
least squares equations of the global fit have to be calculated.
The first contribution is from the first order derivatives with respect
to the global parameters:

\anchor eq-c1 (15)
\f{equation*}{  \label{eq:c1}
\left(\vec{\Delta C}_1\right)_{jk} =
\sum_i \left( \frac{\partial f_i}{\partial p_j} \right)
        \left( \frac{\partial f_i}{\partial p_k} \right)
            \frac{1}{\sigma_i^2}
            \quad \quad \quad \quad
    \Delta  g_j=  \sum_i \left( \frac{\partial f_i}{\partial p_j} \right)
\,  \frac{z_i^{''}}{\sigma_i^2}  \; .
\f}
Note that in the \f$g_j\f$-term the residual is already corrected for the
local fit result. The vector \f$\vec{g}\f$ is the gradient of the global
objective function.

The second contribution is, according to the **Millepede** principle,
a mixed contribution from first order global and local derivatives.
After the local fit the matrix \f$\mathbf{G}\f$ is formed, which has a number of
columns equal to the number of local parameters, and a number of rows
equal to the number of global parameters.
The elements are
\f{equation*}{
\left(\mathbf{G}\right)_{jk} =
\sum_i \left( \frac{\partial f_i}{\partial p_j} \right)
        \left( \frac{\partial f_i}{\partial q_k} \right)
            \frac{1}{\sigma_i^2}    \; .
\f}
The second contribution to the global matrix \f$\mathbf{C}\f$ is then

\anchor eq-c2 (16)
\f{equation*}{ \label{eq:c2}
\vec{\Delta C_2} = - \mathbf{G} \vec{\Gamma}^{-1} \mathbf{G}^{\top}
\f}
with the matrices \f$\mathbf{G}\f$ and \f$\vec{\Gamma}\f$ from a local fit.
Because of the dimension of matrices \f$\mathbf{G}\f$ and  \f$\vec{\Delta C}\f$
the calculation would require a large number of operations. However
only a small number of rows of matrix \f$\mathbf{G}\f$ is not equal to zero
and the calculation can be restricted, with the help of pointers,
to the non-zero part and in addition one can use the fact that the
result is a symmetric matrix. With this code, which may appear
somewhat complicated in comparison to the standard matrix multiplication,
the computation time is small. Note that no extra contribution to the
gradient vector \f$\vec{g}\f$ is necessary because the residual \f$z_i^{''}\f$
used above already includes the local fit result.
The sum \f$\mathbf{C}\f$
of the two contributions,  \f$\vec{\Delta C}_1\f$ from equation \ref eq-c1 "(15)", and
\f$\vec{\Delta C}_2\f$ from equation \ref eq-c2 "(16)", summed over all
local fit-objects, is the final matrix of the least squares
normal equations. Element \f$(\mathbf{C})_{jk}\f$ is related
to the global parameters
with indices \f$j\f$ and \f$k\f$.
The matrix \f$\mathbf{C}\f$ corresponds to a simultaneous fit of the global
*and* local parameters of all local-fit objects. The first term
\ref eq-c1 "(15)"
alone corresponds to the fit, when the local parameters are assumed
to be fixed; only those elements \f$(\mathbf{C})_{jk}\f$ of the matrix \f$\mathbf{C}\f$ are
non-zero, where the two global parameters with indices \f$j\f$ and \f$k\f$
appear in the same measurement \f$y_i\f$. In contrast, in the second term
\ref eq-c2 "(16)" those elements \f$(\mathbf{C})_{jk}\f$ of the matrix \f$\mathbf{C}\f$ are
non-zero, where the two global parameters with indices \f$j\f$ and \f$k\f$
appear in the same local fit.

\paragraph par-downw Outlier downweighting in local fits

The influence of a single measurement \f$y_i\f$ resp. \f$z_i\f$ is proportional
to the absolute value of the normalised residual
\f$ \zeta_i =  z_i/\sigma_i\f$
in the pure least squares fit of local-fit objects.
Thus measurements with a large deviation will in general distort the
resulting fit. This effect can be avoided by downweighting those
measurements. Measurements with large residuals are included in the
method of M-estimates by assigning to them, instead of the ideal
Gaussian distribution of least squares, a different distribution with
larger tails. Technically this is done by multiplying the standard weight
\f$1/\sigma_i^2\f$ of least squares by another factor, which depends on the
actual value of  \f$ \zeta_i = z_i/\sigma_i\f$.
This method requires an initial fit
without downweighting, followed by iterative downweighting.

In **Millepede II** the local fits are improved, if necessary,
by downweighting *after* the first iteration, with a number of iterations
specified by the user with keyword \c outlierdownweighting
{it numberofiterations}. In iterations 2 and 3 the downweighting is done
with the Huber function and the Cauchy function is used, if more iterations
are requested (see section \ref sssec-outlierdeb for an explanation of the
two functions).

\subsubsection sssec-glofit Global fit

\paragraph par-glonocon Global fit without constraints

The aim is to find a step vector \f$\vec{\Delta p}\f$, which minimizes the
objective function: \f$F(\vec{p} + \vec{\Delta p}) = \f$ minimum.
In the \f$k\f$-th iteration the approximate solution is \f$\vec{p}_k\f$.
A quadrative model function \f$\widetilde{F}(\vec{p} + \vec{d})\f$

\anchor eq-quadapp (17)
\f{equation*}{ \label{eq:quadapp}
\widetilde{F}(\vec{p}_k + \vec{d}  ) = F_k + \vec{g}^{\top} \vec{d} +
        \tfrac{1}{2} \vec{d}^{\top}  \mathbf{C} \vec{d}
\f}
is defined with a vector \f$\vec{d}\f$.
The step vector \f$\vec{d}\f$ is determined as solution of the linear system
\f{equation*}{ \label{eq:stepd}
                \mathbf{C} \vec{d} = - \vec{g} \; .
\f}
The step \f$\vec{d}\f$ will  minimize the quadratic function
\f$\widetilde{F}(\vec{p} + \vec{d})\f$ of equation \ref eq-quadapp "(17)",
unless there is some inaccuracy due to rounding errors and other
approximations of the matrix or the solution method. In order
to get a *sufficient* decrease of the objective function
\f$F(\vec{p} + \vec{\Delta p})\f$ a line search algorithm is used.

\paragraph par-linesearch Line search

In the line search algorithm a search is made for the minimum of the
function \f$F(\vec{p} + \vec{\Delta p})\f$ along a line
\f$\vec{p} + \alpha \cdot \vec{d}\f$ in the parameter space, where \f$\alpha\f$
is a factor to be determined. A function \f$\Phi(\alpha)\f$ of this factor,
\f{equation*}{
    \Phi(\alpha) \equiv  F(\vec{p} + \alpha \cdot \vec{d})
\f}
is considered. Each function evaluation at a certain value of \f$\alpha\f$
requires one loop through the data with all local fits. The function
value for \f$\alpha=0\f$ is already calculated before (start value);
the first value
of \f$\alpha\f$ to be considered is \f$\alpha=1\f$ and this value is often
already a sufficiently accurate minimum. In the line search algorithm
the so-called strong Wolfe conditions are used; these are
\f{align*}{
        F( \vec{p} + \alpha \cdot \vec{d}) & \le F( \vec{p})
                + C_1 \alpha \nabla F^{\top} \vec{d} \\
    \left| \nabla F(\vec{p} + \alpha \cdot \vec{d})^{\top} \vec{d}   \right|
        & \le C_2 \left| \nabla F^{\top} \vec{d} \right|
\f}
with \f$0 < C_1 < C_2 < 1\f$. The first condition requires a sufficient decrease
of the function. The second condition, called curvature condition, ensures
a slope reduction, and rules out unacceptable short steps.
In practice the constant \f$C_1\f$ is chosen to small,
for example  \f$C_1 = 10^{-4}\f$. Typical values of the constant  \f$C_2\f$ are
0.9 for a search vector by a Newton method and 0.1 for a search vector
by the conjugate gradient method. Default values in **Millepede II** are
\f$C_1 = 10^{-4}, \, C_2 = 0.9\f$. Often only one value of \f$\alpha\f$,
sometimes two to three values are evaluated. The last function value
should never be higher than the start value.

\paragraph par-iter Iterations

The global fit can be performed with a result \f$\vec{\Delta p}\f$ in one step,
if the matrix equation is accurately solved and if there are no outliers.
The result \f$\vec{\Delta p}\f$ calculated in the first data loop is usually
already close to the final solution; it gives the largest reduction
of the objective function. Outliers introduce some non-linearity into
the problem and require iterations with repeated data loops.
Another argument for repeating the steps is the non-neglectigible
inaccuracy of the numerical solution of the step calculation due to
rounding errors for the large number of parameters.

<b>Iteration 0.</b> The first data loop is called iteration 0. The
function value, the first derivative vector \f$\vec{g}\f$ and the second
derivative matrix \f$\mathbf{C}\f$ or an approximation to it are calculated.
They allow to calculate the first step \f$\vec{\Delta p}= \vec{d}\f$.

<b>Iteration \f$\ge 1\f$.</b> In all further iterations a line search
is performed. Starting value is the point in parameter space, obtained
in the previous iteration. Each function evaluation requires a data loop
with local fits. The expected decrease \f$\Delta F\f$ in an iteration can be
estimated by
\f{equation*}{
            \Delta F_{\textrm{estimate}} = - \vec{g}^{\top} \vec{d} \; .
\f}
This can be compared with the actual decrease
\f{equation*}{
            \Delta F_{\textrm{actual}} = F(\vec{p}) - F(\vec{p}+ \alpha \vec{d})
\f}
at the end of the iteration. Both  values should become smaller during
the iterations; convergence can be assumed
if both are of the order of \f$\approx 1\f$.

\paragraph par-glowithcon Global fit with constraints

In the Lagrange multiplier method  one additional parameter \f$\lambda\f$ is
introduced for each single constraint, resulting in an \f$m\f$-vector
\f$\vec{\lambda}\f$ of Lagrange multiplier. A term depending on \f$\vec{\lambda}\f$
and the constraints is added to the objective function, resulting
in the Lagrange function
\f{equation*}{
\mathcal{L}(\vec{p},\vec{\lambda}) =
F_k + \vec{g}^{\top} \vec{d} +
        \tfrac{1}{2} \vec{d}^{\top}  \mathbf{C} \vec{d}
+ \vec{\lambda}
\left(  \mathbf{A} \left( \vec{p}_k + \vec{d}\right) - \vec{c} \right)
\f}
Using as before a quadratic model for the function \f$F(\vec{p})\f$ and taking
derivatives w.r.t. the parameters \f$\vec{p}\f$ and the Lagrange multipliers
\f$\vec{\lambda}\f$, the two equations
\f{alignat*}{{2}
\mathbf{C} \; & \vec{d} +  \mathbf{A}^{\top} & \vec{\lambda} & = - \vec{g} \\
\mathbf{A} \; &  \vec{d}  &                                 & =
\vec{c} - \mathbf{A} \vec{p}_k
\f}
are obtained; the second of these equations is the constraint equation.
This system of two equations can be combined into one matrix equation
\f{equation*}{ \label{eq:lageq}
\left(
    \begin{array}{ccc|c}
&           & &         \\
&    \mathbf{C} & & \mathbf{A}^{\top}  \\
&           & &         \\   \hline
&     \mathbf{A} & & \mathbf{0}
    \end{array}
    \right)
\left( \begin{array}{c}
~\\   \vec{d}  \\  ~\\  \hline \vec{\lambda}
        \end{array} \right)  =
\left( \begin{array}{c}
~\\  -  \vec{g}  \\  ~\\  \hline   \vec{c} - \mathbf{A} \vec{p}_k
        \end{array} \right)  \; .
\f}
The matrix on the left hand side is still symmetric.
*Linear* least squares problems with *linear*
constraints can be solved directly, without iterations and without
the need for initial values of the parameters. Insufficient
precision of the equality constraints
is corrected by the method discussed in section \ref sssec-feas.

\subsubsection sssec-outfiles Output text files

Several output text file are generated in **Pede**. All files have
a fixed name. If files with this name are existing at the **Pede**
start, they are renamed by appending a \f$\sim\f$ to the file name; existing
files with file name already carrying the \f$\sim\f$ are removed!

<b>Log-file for</b> **Pede**. The file \c millepede.log
summarizes the job execution.
Especially information on the iterations is given.

<b>Result file.</b> The file \c millepede.res contains one line
per (variable or fixed) global parameter; the first three numbers
have the same meaning as the lines following the \c Parameter
keyword in the input text files. The first value is the label, the
second is the fitted parameter value, and the third is the presigma,
either copied from the input text file or zero (default value).

<b>Debug file.</b> The file \c mpdebug.txt, if requested by the
\c printrecord keyword, contains detailed information of all data
and fits for local-fit objects.

<b>Eigenvector file.</b> The file \c millepede.eve is written for the
<tt>method diagonalization</tt>. It contains selected eigenvectors and
their eigenvalues, and can be used e.g. for an analysis of weak modes.

\subsection ssec-align Using Millepede II for track based alignment

\subsubsection sssec-alignpar Global parameters

Global parameters to be determined in a track-based alignment are mostly
geometrical corrections like translation shifts and rotation angles.
Usually only small corrections to the default position and angle information
are necessary, and the assumption of *constant* first derivatives
is sufficient (otherwise an iteration with creation of new data files
would be necessary). In general an alignment should be done simultaneously
of *all* relevant detectors, perhaps using structural constraints
(see below).

In addition to the geometrical corrections there are several additional
parameters, which determine track-hit positions like Lorentz-angle and
\f$T_0\f$-values, drift-velocities etc. These parameters should be included
in the alignment process, since, due to the correlation between all pairs
of global parameters, a
separate determination of those parameters cannot be optimal. Also
beam parameters like vertex position and beam direction angles should be
included in the alignment, if those value are used with the tracks.
Incorrectness of the model used e.g. to calculate the expected
track hit positions, for example a wrong value of the Lorentz-angle,
can introduce distortions in other parameters.

\subsubsection sssec-aligncon Constraints

Equality constraints are, in the mathematical sense, equations
between global parameters, which have to be exact (within the
computing accuracy). Only linear constraints are supported.
Below two areas of constraints are discussed.

\paragraph par-alignundef Removal of undefined degrees of freedom

A general linear transformation has 12 parameters. Three parameters
of the translation and additional three parameters of a rotation
are undefined, if only track residuals are minimized. At least
these six parameters have to be fixed by constraining to zero the overall
translation and rotation.

\paragraph par-alignstruc Structural constraints

A large track detector has usually a sub-detector structure: the
track detector is built from smaller units, which are built by smaller
sub-units etc. down to the level of a single sensor. The different
levels can be visualized by a tree structure. Figure \ref fig-xms "2" shows
the tree structure of the cms detector.

The tree structure can
be reflected in the label structure and constraints.
One set of global parameters is assigned to a unit. To each of the
subunits another set of global parameters; the sum of the transformations
of the subunits is forced to zero by sum constraints, which remove the
otherwise undefined degrees of freedom. For a unit with \f$N\f$ subunits,
with for example 6 degrees of freedom, the total number of degrees of
freedom becomes
\f{equation*}{
\left[ 6 + 6 \cdot N \right] \; \textrm{(parameters)}
- 6 \; \textrm{(constraints)}  = 6 \cdot N  \; .
\f}
The parameters of the unit are more accurately defined than the
parameters of the sub-units. Additional external measurement information
from e.g. a survey can be assigned to the parameters of the unit.
If all sub-unit parameters are fixed after an initial alignment,
the parameters of the unit can still be variable and be improved in
an alignment with a reduced total number of parameters.

\anchor fig-xms Figure 2
\image html structure.png "Figure 2: The tree of components of the cms inner track detectors."
\latexonly
\begin{figure}[h]
\begin{center}
\includegraphics[width=15cm]{img/structure.pdf}
\caption*{Figure 2: The tree of components of the cms inner track detectors.
\label{fig:xms}}
\end{center}
\end{figure}
\endlatexonly

\subsubsection sssec-lintran Linear transformations in 3D

<b>Nominal transformation.</b>
For track reconstruction the local coordinate system and the global
\f$x \, y \,z\f$ coordinate system  are used.
The local system with coordinate \f$\vec{q}\f$
and components \f$u\f$, \f$v\f$ and \f$w\f$ is defined w.r.t. a sensor
with origin at the center of the sensor; the \f$u\f$-axis is along the
precise and the \f$v\f$-axis is along the coarse coordinate, and the
\f$w\f$-axis is normal to the sensor. The transformation from the
global to the local system is given by
\f{equation*}{
\vec{q} = \mathbf{R} \left( \vec{r} - \vec{r}_0 \right)
\f}
with the definitions
\f{equation*}{
\vec{q} =  \left( \begin{array}{c}
u \\ v \\ w \end{array} \right)
\quad \quad \quad
\vec{r} =  \left( \begin{array}{c}
x \\ y \\ z \end{array} \right)
\quad \quad \quad
\vec{r}_0 =  \left( \begin{array}{c}
x_0 \\ y_0 \\ z_0 \end{array} \right)  \; .
\f}
The nominal position \f$\vec{r}_0\f$ and the  rotation matrix \f$\mathbf{R}\f$
are determined by detector assembly and survey information. The convention
is to use the three Euler angles \f$\vartheta\f$, \f$\psi\f$ and \f$\varphi\f$ to
determine the nominal rotation matrix \f$\mathbf{R}\f$:
\f{equation*}{
\mathbf{R} =
\left( \begin{array}{ccc}
\cos \psi \, \cos \varphi - \cos \vartheta \, \sin \psi \, \sin \varphi &
- \cos \psi \, \sin \varphi - \cos \vartheta \, \sin \psi \, \cos \varphi &
\sin \vartheta \, \sin \psi \\
\sin \psi \, \cos \varphi  + \cos \vartheta \, \cos \psi \,  \sin \varphi &
- \sin \psi \,  \sin \varphi + \cos \vartheta \,  \cos \psi \, \cos \varphi &
- \sin \vartheta \, \cos \psi \\
\sin \vartheta \, \sin \varphi &
\sin \vartheta \, \cos \varphi &
\cos \vartheta
\end{array} \right)
\f}
or
\f{equation*}{
\mathbf{R} =
\left( \begin{array}{ccc}
\cos \psi \, \cos \varphi - \cos \vartheta \, \sin \psi \, \sin \varphi &
\sin \psi \, \cos \varphi  + \cos \vartheta \,
\cos \psi \,  \sin \varphi &
\sin \vartheta \, \sin \varphi \\
- \cos \psi \, \sin \varphi - \cos \vartheta \, \sin \psi \, \cos \varphi &
- \sin \psi \,  \sin \varphi + \cos \vartheta \,
\cos \psi \, \cos \varphi &
\sin \vartheta \, \cos \varphi \\
\sin \vartheta \, \sin \psi &
- \sin \vartheta \, \cos \psi &
\cos \vartheta
\end{array} \right)
\f}
The ranges of the angles are:\f$0 \le \vartheta < \pi\f$,
\f$0 \le \psi < 2 \pi\f$, and \f$0 \le \phi < 2 \pi\f$.

A different convention is to parametrize the rotation by
Euler angles \f$\varphi\f$, \f$\vartheta\f$ and \f$\psi\f$ with the rotation
\f$\mathbf{R} = \mathbf{R}_3(\varphi)  \mathbf{R}_2(\vartheta)  \mathbf{R}_3(\psi)\f$,
where \f$R_i(\delta)\f$ is a rotation by an angle \f$\delta\f$ about the axis
\f$\vec{n}_i\f$. The rotation is
\f{equation*}{
\mathbf{R} =
\left( \begin{array}{ccc}
\cos \varphi \cos \vartheta \cos \psi - \sin \varphi \sin \psi &
-\sin \varphi \cos \psi - \cos \varphi \cos \vartheta \sin \psi &
\cos \varphi \sin \vartheta \\
\cos \varphi \sin  \psi + \sin \varphi \cos \vartheta \cos \psi &
\cos \varphi \cos  \psi - \sin \varphi \cos  \vartheta \sin \psi &
\sin \varphi \sin \vartheta \\
-\sin\vartheta \cos \psi &
\sin \vartheta \sin \psi &
\cos \vartheta
\end{array} \right)
\f}
The ranges of the angles are: \f$0 \le \phi < 2 \pi\f$,
\f$0 \le \vartheta < \pi\f$ and \f$0 \le \psi < 2 \pi\f$.

<b>Alignment correction to the transformation.</b>
The alignment procedure determines a correction to the nominal transformation
by an incremental rotation \f$\Delta \mathbf{R}\f$ and a translation \f$\Delta \vec{r}_0\f$.
The combined translation and rotation becomes
\f{align*}{
        \vec{r}_0 &\to \vec{r}_0 + \Delta \vec{r} \\
        \mathbf{R}   &\to \Delta \mathbf{R} \mathbf{R}  \; .
\f}
The correction matrix \f$\Delta \mathbf{R}\f$ is given by small rotations by
\f$\Delta \alpha\f$, \,\f$\Delta \beta\f$ and \f$\Delta \gamma\f$ around the \f$u\f$-axis, the
(new)  \f$v\f$-axis and the (new)  \f$w\f$-axis:
\f{equation*}{
\Delta  \mathbf{R} =  \mathbf{R}_{\alpha}  \mathbf{R}_{\beta}  \mathbf{R}_{\gamma}
\approx \left(
\begin{array}{ccc}
1 & \Delta \gamma &  - \Delta \beta \\
- \Delta \gamma &       1    &  \Delta \alpha \\
    \Delta \beta & - \Delta\alpha & 1
\end{array}  \right) \; ,
\f}
where the approximation has sufficient accuracy for small angles
\f$\Delta \alpha\f$, \,\f$\Delta \beta\f$ and \f$\Delta \gamma\f$. The position correction
\f$\Delta \vec{r}\f$ transforms to the local system as
\f{equation*}{
\Delta \vec{q} = \Delta \mathbf{R} \, \mathbf{R} \, \Delta \vec{r} \quad \quad \quad
\textrm{with} \quad
\Delta \vec{q} =
\left( \begin{array}{c}
\Delta u \\ \Delta v \\ \Delta w \end{array} \right)
\quad \quad \quad
\Delta \vec{r} =
\left( \begin{array}{c}
\Delta x \\ \Delta y \\ \Delta z\end{array} \right)
\f}
These corrections define the corrected (aligned) transformation
from global to local coordinates:
\f{equation*}{
    \vec{q}^{\textrm{aligned}} = \Delta \mathbf{R} \, \mathbf{R}
            \left( \vec{r} - \vec{r}_0 \right)
- \Delta \vec{q}  \; .
\f}
The task of the alignment by tracks is to determine the correction
angles \f$\Delta \alpha\f$, \f$\Delta \beta\f$ and \f$\Delta \gamma\f$ (and thus the
rotation matrix \f$\Delta \mathbf{R}\f$) and the translation vector \f$\Delta \vec{r}\f$
or vector \f$\Delta \vec{q}\f$ for each individual detector element.

\subsection ssec-ackn Acknowledgement

Markus Stoye has worked with different versions of **Millepede II**
during the development phase in the last two years in
alignment studies for the
inner tracker of the cms experiment.
I would like to thank him for the close collaboration, which was
essential for the development of **Millepede II**.
I would like
to thank also Gero Flucke, who contributed the C\f$^{++}\f$ code, and
Claus Kleinwort for the **Millepede I/II** comparison using
data from the H1 experiment.


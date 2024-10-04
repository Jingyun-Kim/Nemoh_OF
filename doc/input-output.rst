############
Input/Output
############

*********************************
Summary of input and output files
*********************************

Following is the list of the user’s input files and the output files for program, where ``[...]`` denote an optional file and ``<...>`` denotes a variable name:

.. list-table:: Nemoh programs input and output files
   :header-rows: 1
   :name: tab:input_output

   *  - Program
      - Input
      - Output
   *  - ``mesh``
      -  * ``Mesh.cal``
         * ``<geomInput-file>``
      -  * ``<mesh-file>``
         * ``mesh/KH.dat``
         * ``mesh/GC_hull.dat``
         * ``mesh/Hydrostatics.dat``
         * ``mesh/Inertia_hull.dat``
         * ``mesh/Description_Full.tec``
         * ``mesh/Description_Wetted.tec``
   *  - ``preProc``
      -  * ``Nemoh.cal``
         * ``<mesh-file>``
      -  * ``Normalvelocities.dat``
         * ``results/FKForce.dat``
         * ``results/FKForce.tec``
         * ``results/index.dat``
         * ``mesh/Mesh.tec``
         * ``mesh/L10.dat``
         * ``mesh/L12.dat``
         * ``mesh/Integration.dat``
         * ``mesh/Kochin.dat``
         * ``mesh/Freesurface.dat``
   *  - ``hydrosCal``
      -  * ``Nemoh.cal``
         * ``Mesh.cal``
         * ``mesh/L10.dat``
         * ``mesh/L12.dat``
      -  * ``mesh/KH.dat``
         * ``mesh/GC_hull.dat``
         * ``mesh/Inertia_hull.dat``
         * ``mesh/Hydrostatics.dat``
         * ``Mechanics/Inertia.dat``
         * ``Mechanics/Kh.dat``
         * [``Mechanics/Km.dat``]
         * [``Mechanics/Badd.dat``]
   *  -  ``solver``
      -  * ``Nemoh.cal``
         * ``input_solver.txt``
         * ``mesh/L10.dat``
         * ``mesh/L12.dat``
      -  * ``logfile.txt``
         * ``results/Forces.dat``
         * [``results/pressure.<id>.dat``]
         * [``results/Kochin.<id>.dat``]
         * [``results/freesurface.<id>.dat``]
         * [``results/sources.<id>.dat``]
   *  - ``postProc``
      -  * ``Nemoh.cal``
         * ``results/Forces.dat``
         * ``results/index.dat``
         * ``results/FKForce.tec``
         * [``Mechanics/Inertia.dat``]
         * [``Mechanics/Kh.dat``]
         * [``Mechanics/Km.dat``]
         * [``Mechanics/Badd.dat``]
      -  * ``results/CA.dat``
         * ``results/CM.dat``
         * ``results/Fe.dat``
         * ``results/ExcitationForce.tec``
         * ``results/DiffractionForce.tec``
         * ``results/RadiationCoefficients.tec``
         * [``results/IRF.tec``]
         * [``results/IRF_excForce.tec``]
         * [``Motion/RAO.dat``]
   *  - ``QTFpreProc``
      -  * ``Nemoh.cal``
         * ``mesh/L10.dat``
         * ``mesh/L12.dat``
         * [``<FS-mesh-file>``]
         * ``Mechanics/Inertia.dat``
         * ``Mechanics/Kh.dat``
         * ``Mechanics/Km.dat``
         * ``Mechanics/Badd.dat``
         * ``results/sources.<id>.dat``
      -  ``QTFPreprocOut/*.bin``
   *  - ``QTFsolver``
      -  * ``Nemoh.cal``
         * ``mesh/L10.dat``
         * ``mesh/L12.dat``
         * ``results/sources/sources.<id>.dat``
         * ``QTFPreprocOut/*.bin``
      -  ``results/QTF/*.dat``
   *  - ``QTFpostproc``
      -  * ``Nemoh.cal``
         * ``results/QTF/*.dat``
      -  * ``results/QTF/OUT_QTFM_N.dat``
         * ``results/QTF/OUT_QTFP_N.dat``

The files mentioned in the previous table are used as follows:

- **User input files:**
   * ``Nemoh.cal``: contains all the computation parameters
   * ``Mesh.cal``: contains geometry information for the ``mesh`` program and the CoG position for ``hydrosCal``
   * ``<mesh-file>``: main mesh of the body
   * ``<geomInput-file>``: geometry of the body for the ``mesh`` program
   * ``Mechanics/Inertia.dat``: rigid body inertia, approximated by Nemoh if not provided
   * ``Mechanics/Kh.dat``: hydrostatic stiffness matrix, computed by Nemoh if not provided
   * ``Mechanics/Km.dat``: added linear stiffness matrix (e.g. to account for mooring), set to zero by Nemoh if not provided
   * ``Mechanics/Badd.dat``: added linear damping matrix (e.g. to account for viscous damping), set to zero by Nemoh if not provided
   * ``input_solver.txt``: contains numerical parameters for the solver
   * ``<FS-mesh-file>``: free-surface mesh file for computation of the FS contribution term in the QTF
- Intermediate files produced by Nemoh:
   * ``Normalvelocities.dat``
   * ``mesh/Integration.dat``: storage of :math:`\vec{n}.dS` and :math:`\vec{x}\wedge\vec{n}.dS` of the integration of forces on the body
   * ``mesh/Kochin.dat``: list of directions for Kochin function computations
   * ``mesh/Freesurface.dat``: Free-surface elevation output mesh
   * ``mesh/L10.dat``: mesh (basically a copy of ``<mesh-file>``)
   * ``mesh/L12.dat``: additional mesh information: center, normal vector and area of panels
   * ``results/index.dat``: summary of computation parameters (degrees of freedom, wave frequencies and directions...)
   * ``results/sources/sources.<id>.dat``: sources intensity on the body for problem ``<id>``, used for QTF
   * ``QTFPreprocOut/*.bin``
- Output files:
   * ``logfile.txt``: solver log file
   * ``mesh/Mesh.tec``: mesh as understood by Nemoh, with panel centers, normal vectors and areas, Tecplot format
   * ``mesh/KH.dat``: hydrostatic stiffness matrix
   * ``mesh/Inertia_hull.dat``: inertia matrix of the hull (assuming the mass is distributed evenly on the shell)
   * ``mesh/Hydrostatics.dat``: summary of hydrostatics computation (buoyancy center, volume displacement, waterplane area, mass)
   * ``mesh/GC_hull.dat``: center of mass of the hull (assuming the mass is distributed evenly on the shell)
   * ``results/FKForce.dat``: Froude-Kylov forces
   * ``results/Fe.dat``: excitation forces, Tecplot format
   * ``results/CM.dat``: added mass coefficients
   * ``results/CA.dat``: radiation damping coefficients
   * ``results/FKForce.tec``: froude-krylov forces, Tecplot format
   * ``results/ExcitationForce.tec``: excitation forces, Tecplot format
   * ``results/DiffractionForce.tec``: diffraction forces, Tecplot format
   * ``results/RadiationCoefficients.tec``: added mass and radiation damping coefficients, Tecplot format
   * ``results/IRF.tec``: Impulse Response Function for the radiation force, Tecplot format
   * ``results/IRF_excForce.tec``: Impulse Response Function for the excitation force, Tecplot format
   * ``Motion/RAO.dat``: Response Amplitude Operator file
   * ``results/pressure.<id>.dat``: pressure on the body for problem ``<id>``, Tecplot format
   * ``results/freesurface.<id>.dat``: free-surface elevation for problem ``<id>``, Tecplot format
   * ``results/Kochin.<id>.dat``: Kochin functions for problem ``<id>``, Tecplot format
   * ``logfileQTF.txt``: QTF pre-processor and solver log file
   * ``results/QTF/*.dat``: QTF output files
   * ``results/QTF/OUT_QTFM_N.dat``: total difference-frequency QTF
   * ``results/QTF/OUT_QTFP_N.dat``: total sum-frequency QTF

Detail descriptions of the input/output file formats are developed in the next subsections.

****************
User input files
****************

Main input file ``Nemoh.cal``
=============================

``Nemoh.cal`` contains all the computation parameters with the format detailed in :numref:`tab:NemohCal`.

.. table:: ``Nemoh.cal`` input file
   :name: tab:NemohCal

   =============== ===== ===== ===== ===== ===== ===== ==================================================================
   File contents example                               Signification
   =================================================== ==================================================================
   \--- Environment ---------------------------------- *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   1025\.                                              Fluid density :math:`\rho` :math:`[kg/m^3]`
   9.81                                                Gravitional acceleration :math:`g` :math:`[m/s^2]`
   200\.                                               Water depth :math:`d` :math:`[m]` (``0.`` for infinite depth)
   0\.             0\.                                 Wave measurement location :math:`\mathbf{x}_{eff} = (x_{eff},y_{eff})` :math:`[m]` so that :math:`\Phi_I(\mathbf{x}, \omega)= - i \frac{g}{\omega} f_0(z) \mathrm{e}^{i \mathbf{k} \cdot (\mathbf{x} - \mathbf{x}_{eff})}`
   \--- Description of floating bodies --------------- *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   1                                                   Number of bodies
   \--- Body 1 --------------------------------------- *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   mesh.dat                                            Name of mesh file
   657             610                                 Number of nodes and number of panels in mesh
   6                                                   Number of degrees of freedom
   1               1\.   0\.   0\.   0\.   0\.   0\.   Surge
   1               0\.   1\.   0\.   0\.   0\.   0\.   Sway
   1               0\.   0\.   1\.   0\.   0\.   0\.   Heave
   2               1\.   0\.   0\.   0\.   0\.   -5\.  Roll about a point (here :math:`(0,0,-5.)`)
   2               0\.   1\.   0\.   0\.   0\.   -5\.  Pitch about a point (here :math:`(0,0,-5.)`)
   2               0\.   0\.   1\.   0\.   0\.   -5\.  Yaw about a point (here :math:`(0,0,-5.)`)
   ...                                                 *This line is repeated for each degree of freedom*
   6                                                   Number of resulting generalised forces
   1               1\.   0\.   0\.   0\.   0\.   0\.   Force in x direction
   1               0\.   1\.   0\.   0\.   0\.   0\.   Force in y direction
   1               0\.   0\.   1\.   0\.   0\.   0\.   Force in z direction
   2               1\.   0\.   0\.   0\.   0\.   -5\.  Moment force in x direction about a point (here :math:`(0,0,-5.)`)
   2               0\.   1\.   0\.   0\.   0\.   -5\.  Moment force in y direction about a point (here :math:`(0,0,-5.)`)
   2               0\.   0\.   1\.   0\.   0\.   -5\.  Moment force in z direction about a point (here :math:`(0,0,-5.)`)
   ...                                                 *This line is repeated for each generalised force*
   0                                                   *Reserved, single line*
   ...                                                 *This section (header included) is repeated for each body*
   \--- Load cases ----------------------------------- *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   1               100   0.062 6.28                    Wave frequency range: unit (1 for *rad/s*, 2 for *Hz* and 3 for *s*), :math:`N_{\omega}`, :math:`\omega_{min}>0` and :math:`\omega_{max}`
   2               0     30                            Wave direction range: :math:`N_{\beta}`, :math:`\beta_{min}` and :math:`\beta_{max}` (in :math:`[deg]`)
   \--- Post processing ------------------------------ *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   0               0.1   10\.                          IRF (Impulse Response Function) flag (0/1), time step and duration
   0                                                   Pressure output flag (0/1)
   0               0\.   180\.                         Kochin functions direction range: :math:`{N_{\beta}}_{Kochin}`, :math:`{\beta_{min}}_{Kochin}` and :math:`{\beta_{max}}_{Kochin}` (in :math:`[deg]`)
   0               50    400\. 400\.                   Free surface elevation output: number of points in x direction (0 to deactivate) and y direction and (x,y) dimensions of domain
   1                                                   RAO (Response Amplitude Operator) flag (0/1)
   1                                                   Output requency unit, 1 for *rad/s*, 2 for *Hz* and 3 for *s*
   \-- QTF ------------------------------------------- *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   1                                                   QTF (Quadratic Transfer Function) flag (0/1)
   65              0.062 3.14                          :math:`{N_{\omega}}_{QTF}`, :math:`{\omega_{min}}_{QTF}>0` and :math:`{\omega_{max}}_{QTF}` (in :math:`[rad/s]`)
   1                                                   Bidirectional QTF computation flag (0/1)
   2                                                   Contributing terms: 1 ``DUOK``, 2 ``DUOK+HASBO``, 3 Full QTF (``DUOK+HASBO+HASFS+ASYMP``)
   NA                                                  Name of free surface meshfile (**only for full QTF**), 'NA' if not applicable
   0               0     0                             Free surface QTF parameters: :math:`R_e`, :math:`N_{R_e}` and :math:`N_{Bessel}` (**only for full QTF**)
   0                                                   Include hydrostatic terms of the quadratic first order motion (:math:`-[\boldsymbol K] \tilde{\boldsymbol\xi}^{(2)}`), flag (0/1)
   1                                                   For QTFposProc: output frequency unit, 1 for *rad/s*, 2 for *Hz* and 3 for *s*
   1                                                   For QTFposProc: include ``DUOK`` in total QTF, flag (0/1)
   1                                                   For QTFposProc: include ``HASBO`` in total QTF, flag (0/1)
   0                                                   For QTFposProc (**only for full QTF**): include ``HASFS+ASYMP`` in total QTF, flag (0/1)
   =============== ===== ===== ===== ===== ===== ===== ==================================================================

.. warning::
   It is strongly recommended to use the same point for all rotations and moments when computing the RAO (and QTF).
   The inertia, added damping and added stiffness matrices must also be provided at the same reference point.

.. warning::
   The mathematical problem is not defined for :math:`\omega=0`, hence the lower frequency boundary :math:`\omega_{min}` **must be strictly positive**.
   The same goes for :math:`{\omega_{min}}_{QTF}`.

.. note::
   When the frequencies needed for the QTF were not computed by the 1st order solver (Nemoh1), the values are interpolated.


.. important::
   In the case of QTF computation, it is strongly recommended that the computed first-order hydrodynamic coefficients in Nemoh1 cover all difference-frequencies and sum-frequencies intervals for the QTF computation.

   The suggested radial frequency interval is :math:`\omega \in [\Delta \omega, \omega_{max}]` with a step :math:`\Delta \omega` and :math:`\omega_{max}=N_{\omega}\Delta \omega`, where :math:`N_{\omega}` is the number of radial frequencies.

   :math:`\omega_{max}` should be chosen as the maximum computed sum-frequencies :math:`\omega_1+\omega_2`, that is :math:`{\omega_{max}}_{QTF}\leq\omega_{max}/2`.


The meaning of the QTF contributing terms parameter is as follows:

1. ``DUOK``: only compute the quadratic terms of QTFs

2. ``DUOK+HASBO``: compute the quadratic and the body force contributions in the potential QTFs

3. ``DUOK+HASBO+HASFS+ASYMP``: includes the computation of the free-surface integrals in the finite and semi-infinite domains.

The computation of the ``HASFS+ASYMP`` terms requires a free-surface mesh, whose name is specified in ``Nemoh.cal``, along with the following parameters:

-  :math:`R_e`: external radius of the free-surface mesh

-  :math:`N_{R_e}`: number of radial discretization points between the waterline and the external radius :math:`R_e`

-  :math:`N_{Bessel}`: the number of Bessel functions (use ``30`` as a default value)

The hydrostatic terms of the quadratic first-order motion :math:`-[\boldsymbol K] \tilde{\boldsymbol\xi}^{(2)}` may be added to the quadratic QTF term (DUOK), where

   .. math::

      \begin{aligned}
      \tilde{\boldsymbol\xi}^{(2)}=[0,0,z_G(\theta_1^{(1)2}+\theta_2^{(1)2})/2,\theta_2^{(1)}\theta_3^{(1)}/2,-\theta_3^{(1)}\theta_1^{(1)}/2,0]^T.
      \end{aligned}

with :math:`z_G` the vertical component of CoG. Note that this term is optional and needed only in ``QTFsolver``. This term is not always included in other software, *e.g.* HydroSTAR :cite:p:`HYDROSTAR` does not included it.

``QTFpostproc`` computes the total QTF by summing all the terms specified in ``Nemoh.cal``, provided the required terms were computed.

Mesh file
=========

The ``<mesh-file>`` contains all the mesh information with a format as shown in :numref:`tab:meshfile`. Lid panels (:math:`z=0`) of the structure may be included in this file to activate the irregular frequencies removal method. This mesh file may be generated by Nemoh's ``mesh`` program or by an external mesh generator.

.. warning::

   Only the immersed part of the body (:math:`z<0`) should be described in the mesh file.

   Optionally, for surface-piercing bodies, a "lid" may be added to close the open internal surface at :math:`z=0` and avoid irregular frequencies.

External mesh generators, *e.g.* the open-source software GMSH :cite:p:`GMSH`, may be used to generate mesh files but they must be adapted to the Nemoh format.
Using the Python tool `meshmagick <https://github.com/LHEEA/meshmagick>`_ is recommended for this (``mar`` format).
Alternatively, a Matlab file for converting GMSH mesh file to the Nemoh format is provided in the dedicated directory (see :ref:`getting-started:Supporting Matlab files`).

.. table:: ``meshfile`` format
   :name: tab:meshfile

   ======= ============= ============= ============= ============================================================================================================
   File contents                   Signification
   ================================================= ============================================================================================================
   2       1                                         First column must be a 2. Second column is 1 for a symmetric (about :math:`xOz`) body half-mesh, 0 otherwise.
   1       :math:`x_1`   :math:`y_1`   :math:`z_1`   Table of nodes: first column is the node ID, other 3 are the coordinates :math:`(x,y,z)` of each node, listed as rows.
   ...     ...           ...           ...
   0       0\.           0\.           0\.           Last line of the table of nodes.
   1       2             3             4             Table of connectivities: node IDs of each panel listed as rows.
   ...     ...           ...           ...
   0       0             0             0             Last line of the table of connectivities.
   ======= ============= ============= ============= ============================================================================================================


Geometry meshing input file ``Mesh.cal``
========================================

``Mesh.cal:`` contains mesh and environmental parameters with a format as in :numref:`tab:meshcal`. This file is used as input for ``mesh`` and ``hydrosCal``.
All the parameters are used in ``mesh``. Only center of gravity, water density, and gravity are used in ``hydrosCal``.

.. table:: ``Mesh.cal`` file format
   :name: tab:meshcal

   =============== === === ==================================================================
   File contents           Signification
   ======================= ==================================================================
   geomInput_name          Name of the geomInput file.
   0                       1 for a symmetric (about :math:`xOz`) body half-mesh, 0 otherwise.
   0\.             0\.     Translation about x and y axis (respectively)
   0\.             0\. -7  Coordinates of gravity centre
   500\.                   Target for the number of panels in refined mesh
   2\.                     Minimum subdivision of a geometric panel
   0\.
   1\.                     Scaling factor
   1025                    Water density :math:`(kg/m^3)`
   9.81                    Gravity acceleration :math:`(m/s^2)`
   =============== === === ==================================================================

Geometry file
=============

``<geomInput-file>`` contains a geometric description of the body in the form of a (very) coarse mesh, that are number of nodes, number of panels, table of nodes and table of connectivities.
The input file has to follow the format as shown in :numref:`tab:geomInput`.
This format can also be converted from other popular mesh formats using `meshmagick <https://github.com/LHEEA/meshmagick>`_ (``nem`` format).

.. note::

   The geometry file may include parts of the body above the free surface, which will be cut out by the ``mesh`` program.

.. table:: ``geomInput`` file format
   :name: tab:geomInput

   ============= ============= ============= ==== ==================================================================
   File contents                Signification
   ============================================== ==================================================================
   100                                            Total number of nodes.
   25                                             Total number of panels.
   :math:`x_1`   :math:`y_1`   :math:`z_1`        Table of nodes: coordinates :math:`(x,y,z)` of each node listed as rows.
   ...           ...           ...           ...
   1             2             3             4    Table of connectivities: node IDs of each panel listed as rows.
   ...           ...           ...           ...
   ============= ============= ============= ==== ==================================================================


Numerical parameters input file ``input_solver.txt``
====================================================

``input_solver.txt`` contains solver parameters with format as in :numref:`tab:input_solver`. The parameters are described as follows.

-  Number of Gauss Quadrature points, :math:`N^2`, is used for the surface integration in the influence coefficients. User specifies an integer value of :math:`N\in [1,4]`, default :math:`N=2`.

-  Minimum z of flow and source points is defined with a factor :math:`\epsilon_{zmin}` multiplied by the maximal horizontal distance between two point of the mesh, default :math:`\epsilon_{zmin}=0.001`.

-  Three linear-system solvers are available; 1 Gauss elimination, 2 LU Decomposition, 3 GMRES iterative solver.

-  If GMRES solver is chosen then the three parameters, the restart parameter, the relative tolerance and the maximum number of iterations, have to be specified. If the tolerance is not achieved after the maximum iteration exceeded then LU decomposition solves the system directly.

.. table:: ``input_solver.txt`` file format
   :name: tab:input_solver

   =========== ===== ===== =====================================================================================
   File contents           Signification
   ======================= =====================================================================================
   2                       Gauss quadrature order N=\[1,4\] for surface integration, resulting in :math:`N^2` nodes.
   0.001                   eps_zmin for determining minimum z of flow and source points of panel.
   1                       Solver option: 0 for GAUSS ELIM., 1 for LU DECOMP., 2 for GMRES.
   10          1e-5  1000  GMRES parameters: restart parameter, relative tolerance and max number of iterations.
   =========== ===== ===== =====================================================================================



Free surface mesh file
======================

``FSmeshfile`` contains all the free-surface mesh information with a format as shown in :numref:`tab:FSmeshfile`. Quadrilateral panels discretized free-surface area in between the body waterline, :math:`R_B`, and the exterior radius :math:`R_e`. Waterline on :math:`R_B` and :math:`R_e` has to discretized by line segments.

.. table:: ``FSmeshfile`` format (Free surface mesh file)
   :name: tab:FSmeshfile

   ======= ============= ============= ============= ============================================================================================================
   File contents                   Signification
   ================================================= ============================================================================================================
   1       5000          4900           400          Free-surface computation parameters: first column is 1 for a symmetric (about :math:`xOz`) body half-mesh, 0 otherwise. Column 2-4 are number of nodes, number of panels and number of segments for the waterline, respectively.
   1       :math:`x_1`   :math:`y_1`   :math:`z_1`   Table of nodes: first column is the node ID, other 3 are the coordinates :math:`(x,y,z)` of each node, listed as rows.
   ...     ...           ...           ...
   0       0\.           0\.           0\.           Last line of the table of nodes.
   1       2             3             4             Table of connectivities: node IDs of each panel listed as rows.
   ...     ...           ...           ...
   4901    4902                                      Table of connectivities for the waterline: node IDs of each segment listed as rows.
   ...     ...           ...           ...
   0       0             0             0             Last line of the table of connectivities.
   ======= ============= ============= ============= ============================================================================================================


Matrix files
============

``Km.dat`` and ``Badd.dat`` are additional stiffness matrix and damping coefficient matrix. The files contains the matrix components with size :math:`(N_{body}\cdot N_{DoF})\times (N_{body}\cdot N_{DoF})`.


************
Output files
************

Inertia and hydrostatic stiffness
=================================

Hydrostatic output files such as inertia and stiffness matrices are produced by ``mesh``, if ``<geomInput-file>`` is prescribed, or by ``hydrosCal``, if ``meshfile`` is prescribed. The files contain the matrix components with size :math:`(N_{body}\cdot N_{DoF})\times (N_{body}\cdot N_{DoF})`.

Forces
======

``FKForce.tec``, ``DiffractionForce.tec`` and ``ExcitationForce.tec`` are the output files of the Froude-Krylov, the diffraction and the excitation forces respectively. The output file format is given in :numref:`tab:WaveForce`.
The file contains the absolute value and the phase [deg] of the force for each ’frequency’ :math:`\omega`. The force is given for each specified force axis (i.e. surge, heave, pitch) for each body.
The ’frequency’ is given based on the chosen type, [rad/s, Hz, s], of the post-processing parameter in ``Nemoh.cal``, except the Froude-Krylov force, which is only in the radial frequency [rad/s].

.. table:: Output file format of Froude-Krylov, diffraction and excitation forces
   :name: tab:WaveForce

   =========================== ================================== ======================================= ================ ================ ============================================= ================================================
   :math:`\omega_1`            :math:`|F_1(\omega_1)|`            :math:`\angle F_1(\omega_1)`            :math:`\cdots`   :math:`\cdots`   :math:`|F_{N_{forces}}(\omega_1)|`            :math:`\angle F_{N_{forces}}(\omega_1)`
   :math:`\omega_2`            :math:`|F_1(\omega_2)|`            :math:`\angle F_1(\omega_2)`            :math:`\cdots`   :math:`\cdots`   :math:`|F_{N_{forces}}(\omega_2)|`            :math:`\angle F_{N_{forces}}(\omega_2)`
   :math:`\vdots`              :math:`\vdots`                     :math:`\vdots`                          :math:`\vdots`   :math:`\vdots`   :math:`\vdots`                                :math:`\vdots`
   :math:`\omega_{N_\omega}`   :math:`|F_1(\omega_{N_\omega})|`   :math:`\angle F_1(\omega_{N_\omega})`   :math:`\cdots`   :math:`\cdots`   :math:`|F_{N_{forces}}(\omega_{N_\omega})|`   :math:`\angle F_{N_{forces}}(\omega_{N_\omega})`
   =========================== ================================== ======================================= ================ ================ ============================================= ================================================


Radiation coefficients
======================

``RadiationCoefficients.tec`` is the output file for added mass and damping coefficients with format as in :numref:`tab:addedmass_damping_coeffs`. The radiation coefficients are given for each :math:`DoF`, each force axis and for each frequency. The frequency is given based on the chosen frequency/period unit (*rad/s*, *Hz* or *s*) in ``Nemoh.cal``.

The hydrodynamic coefficients are also produced in the *.dat* files, i.e. *CA.dat* for the damping coefficients, *CM.dat* for the added mass coefficients, *Fe.dat* for the excitation force and *FKForce.dat* for the excitation force. The frequency type of the output files is only radial frequency [rad/s]. These output files are used as input files for the QTF modulus.

.. table:: Output file format of the radiation coefficients
   :name: tab:addedmass_damping_coeffs

   ============================== =========================================== ======================================== ================ ================ ==================================================== ====================================================
   :math:`\omega_1`               :math:`M^a_{11}(\omega_1)`                  :math:`B_{11}(\omega_1)`                 :math:`\cdots`   :math:`\cdots`   :math:`M^a_{1N_{forces}}(\omega_1)`                  :math:`B_{1N_{forces}}(\omega_1)`
   :math:`\omega_2`               :math:`M^a_{11}(\omega_2)`                  :math:`B_{11}(\omega_2)`                 :math:`\cdots`   :math:`\cdots`   :math:`M^a_{1N_{forces}}(\omega_2)`                  :math:`B_{1N_{forces}}(\omega_2)`
   :math:`\vdots`                 :math:`\vdots`                              :math:`\vdots`                           :math:`\vdots`   :math:`\vdots`   :math:`\vdots`                                       :math:`\vdots`
   :math:`\omega_{N_\omega}`      :math:`M^a_{11}(\omega_{N_\omega})`         :math:`B_{11}(\omega_{N_\omega})`        :math:`\cdots`   :math:`\cdots`   :math:`M^a_{1N_{forces}}(\omega_{N_\omega})`         :math:`B_{1N_{forces}}(\omega_{N_\omega})`
   :math:`\omega_1`               :math:`M^a_{21}(\omega_1)`                  :math:`B_{21}(\omega_1)`                 :math:`\cdots`   :math:`\cdots`   :math:`M^a_{2N_{forces}}(\omega_1)`                  :math:`B_{2N_{forces}}(\omega_1)`
   :math:`\vdots`                 :math:`\vdots`                              :math:`\vdots`                           :math:`\vdots`   :math:`\vdots`   :math:`\vdots`                                       :math:`\vdots`
   :math:`\omega_{N_\omega}`      :math:`M^a_{21}(\omega_{N_\omega})`         :math:`B_{21}(\omega_{N_\omega})`        :math:`\cdots`   :math:`\cdots`   :math:`M^a_{2N_{forces}}(\omega_{N_\omega})`         :math:`B_{2N_{forces}}(\omega_{N_\omega})`
   :math:`\vdots`                 :math:`\vdots`                              :math:`\vdots`                           :math:`\vdots`   :math:`\vdots`   :math:`\vdots`                                       :math:`\vdots`
   :math:`\omega_{N_\omega}`      :math:`M^a_{N_{DoF}1}(\omega_{N_\omega})`   :math:`B_{N_{DoF}1}(\omega_{N_\omega})`  :math:`\cdots`   :math:`\cdots`   :math:`M^a_{N_{DoF}N_{forces}}(\omega_{N_\omega})`   :math:`B_{N_{DoF}N_{forces}}(\omega_{N_\omega})`
   ============================== =========================================== ======================================== ================ ================ ==================================================== ====================================================


Response Amplitude Operator
===========================

``RAO.dat`` is the output file of the response amplitude operator with the file format as in :numref:`tab:RAO`. The output file gives the absolute value and the phase of RAO for each degree of freedom and each frequency. The frequency is given based on the chosen ’frequency’ type, [rad/s, Hz, s], of the post-processing parameter in ``Nemoh.cal``. Only radial frequency output file will be produced in the case of the QTF computed.

.. table:: Output file format of ``RAO.dat``
   :name: tab:RAO

   =========================== ==================================== ================ ==================================== ========================================= ================ =============================================
   :math:`\omega_1`            :math:`|\xi_1(\omega_1)|`            :math:`\cdots`   :math:`|\xi_6(\omega_1)|`            :math:`\angle \xi_1(\omega_1)`            :math:`\cdots`   :math:`\angle \xi_6(\omega_1)`
   :math:`\vdots`              :math:`\vdots`                       :math:`\vdots`   :math:`\vdots`                       :math:`\vdots`                            :math:`\vdots`   :math:`\vdots`
   :math:`\omega_{N_\omega}`   :math:`|\xi_1(\omega_{N_\omega})|`   :math:`\cdots`   :math:`|\xi_6(\omega_{N_\omega})|`   :math:`\angle \xi_1(\omega_{N_\omega})`   :math:`\cdots`   :math:`\angle \xi_6(\omega_{N_\omega})`
   =========================== ==================================== ================ ==================================== ========================================= ================ =============================================


Impulse Response Functions
==========================

``IRF.tec`` and ``IRF_excForce.tec`` are the impulse response functions for the radiation damping and the excitation force, respectively. The radiation damping IRF has the file format as in :numref:`tab:IRF` and the excitation force IRF as in :numref:`tab:IRFExcF`.

.. table:: Output file format of ``IRF.tec``
   :name: tab:IRF

   ================== ================================ ============================= ================ =================== ========================================= ====================================
   :math:`t_1`        :math:`M^a_{11}(\infty)`         :math:`IRF_{11}(t_1)`         :math:`\cdots`   :math:`\cdots`      :math:`M^a_{1N_{forces}}(\infty)`         :math:`IRF_{1N_{forces}}(t_1)`
   :math:`t_2`        :math:`M^a_{11}(\infty)`         :math:`IRF_{11}(t_2)`         :math:`\cdots`   :math:`\cdots`      :math:`M^a_{1N_{forces}}(\infty)`         :math:`IRF_{1N_{forces}}(t_2)`
   :math:`\vdots`     :math:`\vdots`                   :math:`\vdots`                :math:`\vdots`   :math:`\vdots`      :math:`\vdots`                            :math:`\vdots`
   :math:`t_1`        :math:`M^a_{21}(\infty)`         :math:`IRF_{21}(t_1)`         :math:`\cdots`   :math:`\cdots`      :math:`M^a_{2N_{forces}}(\infty)`         :math:`IRF_{2N_{forces}}(t_1)`
   :math:`\vdots`     :math:`\vdots`                   :math:`\vdots`                :math:`\vdots`   :math:`\vdots`      :math:`\vdots`                            :math:`\vdots`
   :math:`t_N`        :math:`M^a_{N_{DoF}1}(\infty)`   :math:`IRF_{N_{DoF}1}(t_N)`   :math:`\cdots`   :math:`\cdots`      :math:`M^a_{N_{DoF}N_{forces}}(\infty)`   :math:`IRF_{N_{DoF}N_{forces}}(t_N)`
   ================== ================================ ============================= ================ =================== ========================================= ====================================


.. table:: Output file format of ``IRF_excForce.tec``
   :name: tab:IRFExcF

   ================ ====================== ================ =============================
   :math:`t_1`      :math:`IRF_{1}(t_1)`   :math:`\cdots`   :math:`IRF_{N_{forces}}(t_1)`
   :math:`t_2`      :math:`IRF_{1}(t_2)`   :math:`\cdots`   :math:`IRF_{N_{forces}}(t_2)`
   :math:`\vdots`   :math:`\vdots`         :math:`\vdots`   :math:`\vdots`
   :math:`t_N`      :math:`IRF_{1}(t_N)`   :math:`\cdots`   :math:`IRF_{N_{forces}}(t_N)`
   ================ ====================== ================ =============================


Field results
=============

``pressure.<id>.dat``, ``kochin.<id>.dat`` and ``freesurface.<id>.dat`` are output files for the pressure, Kochin function and free surface respectively, for a specific problem number ``<id>``. The problem number is defined in order of the :math:`N_{\beta}` diffraction problems and the :math:`N_{DoF}` radiation problems, all for each frequency. So problem number ``001`` is the diffraction problem for the first frequency and first wave direction.

The `<id>` of the diffraction problem for frequency number :math:`i_{\omega}` and wave direction number :math:`i_{\beta}` is:

.. math::

   <id> = [N_{\beta} + N_{DoF}]*(i_{\omega} - 1) + i_{\beta}

The `<id>` of the radiation problem for frequency number :math:`i_{\omega}` and degree of freedom number :math:`i_{DoF}` is:

.. math::

   <id> = N_{\beta}i_{\omega} + N_{DoF}(i_{\omega} - 1) + i_{DoF}


For example, if :math:`N_{\beta}=1`, then problem ``002`` is the radiation problem for the first degree of freedom and the first frequency. Additionally if :math:`N_{DoF}=6` then problem ``008`` is the diffraction problem for the second frequency.

Pressure
^^^^^^^^

``pressure.<id>.dat`` is a pressure output file for problem number ``<id>``.
In each file, the absolute value of pressure :math:`|P|` (*Pa*) and the phase :math:`\angle P` (*rad*) are given for each body panel.
The format of the output file is given in :numref:`tab:pressure`.

.. table:: Output file format of ``pressure.<id>.dat``
   :name: tab:pressure

   ======================= ======================= ======================= ======================================== ===========================================
   :math:`x_1`             :math:`y_1`             :math:`z_1`             :math:`|P(\boldsymbol x_1)|`             :math:`\angle P(\boldsymbol x_1)`
   :math:`\vdots`          :math:`\vdots`          :math:`\vdots`          :math:`\vdots`                           :math:`\vdots`
   :math:`x_{N_{panel}}`   :math:`y_{N_{panel}}`   :math:`z_{N_{panel}}`   :math:`|P(\boldsymbol x_{N_{panel}})|`   :math:`\angle P(\boldsymbol x_{N_{panel}})`
   ======================= ======================= ======================= ======================================== ===========================================


Kochin functions
^^^^^^^^^^^^^^^^

``kochin.<id>.dat`` is an output file of the Kochin function on a prescribed direction for for problem number ``<id>``.
In each file, depending on the diffraction/radiation problem, the computed absolute value of the Kochin function :math:`|\mathcal{H}|` and the phase :math:`\angle \mathcal{H}` (*rad*) are saved for each direction :math:`\vartheta`.
The format of the output file is given in :numref:`tab:kochin`.

.. table:: Output file format of *kochin.<id>.dat*
   :name: tab:kochin

   ================================ =============================================== ==================================================
   :math:`\vartheta_1`              :math:`|\mathcal{H}(\vartheta_1)|`              :math:`\angle \mathcal{H}(\vartheta_1)`
   :math:`\vdots`                   :math:`\vdots`                                  :math:`\vdots`
   :math:`\vartheta_{N\vartheta}`   :math:`|\mathcal{H}(\vartheta_{N\vartheta})|`   :math:`\angle \mathcal{H}(\vartheta_{N\vartheta})`
   ================================ =============================================== ==================================================


Free surface elevation
^^^^^^^^^^^^^^^^^^^^^^

``freesurface.<id>.dat`` is an output file of the free-surface elevation on a prescribed free-surface domain for problem number ``<id>``.
In each file, depending on the diffraction/radiation problem, the computed absolute value of the free-surface elevation :math:`|\eta|` and the phase :math:`\angle \eta` (*rad*) are saved for each free-surface position.
The format of the output file is given in :numref:`tab:freesurface`.

.. table:: Output file format of ``freesurface.<id>.dat``
   :name: tab:freesurface

   ======================= ======================= ===================================== ========================================== ======================================== ======================================
   :math:`x_1`             :math:`y_1`             :math:`|\eta(\vec{x}_1)|`             :math:`\angle \eta(\vec{x}_1)`             :math:`Re[ \eta(\vec{x}_1)]`             :math:`Im[ \eta(\vec{x}_1)]`
   :math:`\vdots`          :math:`\vdots`          :math:`\vdots`                        :math:`\vdots`                             :math:`\vdots`                           :math:`\vdots`
   :math:`x_{N_{panel}}`   :math:`y_{N_{panel}}`   :math:`|\eta(\vec{x}_{N_{panel}})|`   :math:`\angle \eta(\vec{x}_{N_{panel}})`   :math:`Re[ \eta(\vec{x}_{N_{panel}})]`   :math:`Im[ \eta(\vec{x}_{N_{panel}})]`
   ======================= ======================= ===================================== ========================================== ======================================== ======================================


Quadratic Transfer Function
===========================

``OUT_QTFM_N.dat`` and ``OUT_QTFP_N.dat`` are the output files of difference- and sum-frequencies QTF.
The QTF results are either the total QTF or parts of the QTF terms that depend on the user choice QTF post-processing parameters in ``Nemoh.cal``.
The QTF values are normalized by :math:`\rho g` and given both as modulus and argument (in *deg*), and as real and imaginary parts.
The frequency/period unit (*rad/s*, *Hz* or *s*) depends on the setting in ``Nemoh.cal``.
The format of the output file is given in :numref:`tab:QTF`. Only the lower triangular part of the QTF matrix is saved in the file.
The upper triangular part can be reconstructed knowing that the QTF matrix is conjugate-symmetric for the difference-frequency and symmetric for the sum-frequency.
A Matlab file for reading this output file is provided in ``matlabRoutines/``.

.. table:: Output file format of ``OUT_QTFM_N.dat`` and ``OUT_QTFP_N.dat``
   :name: tab:QTF

   =============================== ================================= =========================== =========================== ======================= ====================== ==================== ======================== ======================
   :math:`\omega_{1_1}`            :math:`\omega_{2_1}`              :math:`\beta_{1_1}`         :math:`\beta_{2_1}`         :math:`DoF_1`           :math:`|QTF|/\rho g`   :math:`\angle QTF`   :math:`Re[QTF]/\rho g`   :math:`Im[QTF]/\rho g`
   :math:`\omega_{1_2}`            :math:`\omega_{2_1}`              :math:`\beta_{1_1}`         :math:`\beta_{2_1}`         :math:`DoF_1`           :math:`|QTF|/\rho g`   :math:`\angle QTF`   :math:`Re[QTF]/\rho g`   :math:`Im[QTF]/\rho g`
   :math:`\vdots`                  :math:`\vdots`                    :math:`\vdots`              :math:`\vdots`              :math:`\vdots`          :math:`\vdots`         :math:`\vdots`       :math:`\vdots`           :math:`\vdots`
   :math:`\omega_{1_{N_\omega}}`   :math:`\omega_{2_{N_\omega}}`     :math:`\beta_{1_{Nbeta}}`   :math:`\beta_{2_{Nbeta}}`   :math:`DoF_{N_{DoF}}`   :math:`|QTF|/\rho g`   :math:`\angle QTF`   :math:`Re[QTF]/\rho g`   :math:`Im[QTF]/\rho g`
   =============================== ================================= =========================== =========================== ======================= ====================== ==================== ======================== ======================

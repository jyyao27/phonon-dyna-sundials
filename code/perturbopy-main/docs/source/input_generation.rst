Automatic Generation of Input Files for Perturbo
================================================

The PERTURBO code has many calculation modes (specified by the  `calc_mode <https://perturbo-code.github.io/mydoc_param_perturbo.html#calc_mode>`_ variable). Each calculation mode implies different mandatory and optional input parameters. In order to simplify and systematize the input files for the user, Perturbopy package conatin command ``input_generation``, which generates the PERTURBO input files for different calculation modes.

To use the script, just go to the folder, where you want to generate the file, and run the command (suppose, that we generate the input file for the calculation mode `ephmat`):

.. code-block:: console

   $ input_generation --calc_mode ephmat

For a shorter version, one can specify ``-c`` instead of ``--calc_mode``.

Then, the input file (called by default *pert.in*) is generated:

.. code-block:: console

	! This input file for PERTURBO was generated by Perturbopy 
	! Date: November 20, 2024  14:42

	&perturbo
	! ***Mandatory parameters***
	calc_mode = 'ephmat'
	prefix = 'prefix'
	fklist = 'prefix.kpt'
	fqlist = 'prefix.qpt'


	! ***Optional parameters***
	! band_min = 1
	! band_max = 9999999
	! phfreq_cutoff = 1.0
	/
   

.. code-block:: python
	! This input file for PERTURBO was generated by Perturbopy 
	! Date: November 20, 2024  14:42 

	&perturbo
	! ***Mandatory parameters***
	calc_mode = 'ephmat'
	prefix = 'prefix'
	fklist = 'prefix.kpt'
	fqlist = 'prefix.qpt'


	! ***Optional parameters***
	! band_min = 1
	! band_max = 9999999
	! phfreq_cutoff = 1.0
	/

It contains a block of mandatory parameters for this calculation mode and a block of optional ones, which is commented. As one can see, this input file containes some *typical* values for the input parameters. The user should modify them for a given calculation. 

Setting the variables is also possible using the scipt. For example, to set ``prefix`` to `'si'` and ``band_min`` to `'10'`, run:

.. code-block:: console

   $ input_generation -c ephmat --prefix si --band_min 10

The values of these parameters were changed in the *pert.in* file. Note, that since we specified an optional parameter ``band_min``, it was uncommented in the input file.

To change the name of the input file, run the script with ``-i your_name.in`` option. Setting the input parameter values from the scipt could be usefull in the case, when one needs to create automatically many different input files. 

In order to generate the input files for the ``qe2pert.x`` calcuation, select ``-c qe2pert``. Run the script with ``-h`` to get the whole list of possible options.

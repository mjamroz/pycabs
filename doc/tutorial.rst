
=========
Tutorial
=========


Calculating heat capacity, :math:`C_v` 
=======================================
.. math:: 
	C_v(T) = \frac{\left<E^2\right> - \left<E\right>^2}{T^2}


.. literalinclude:: ../examples/heat_capacity.py

.. image:: images/cv.png


Calculating dependency of average number of native contacts on temperature :math:`T`
====================================================================================

Contact in script below is defined as 1 if distance between C-alpha atoms is below :math:`8\AA` and sequence distance > 5 residues. Script run simulations in parallel fashion, and display plot at the end.

.. literalinclude:: ../examples/transition_temp.py

.. image:: images/avg_nat_cont.png


Monitoring of CABS energy during simulation
===========================================

.. literalinclude:: ../examples/monitoring_energy.py

Monitoring of end-to-end distance of chain during simulation
============================================================

.. literalinclude:: ../examples/monitoring_e2e_distance.py

Detector Servers
=================

Location
---------
   slsDetectorPackage/serverBin/ folder in every release.


.. _Detector Server Arguments:

Arguments
---------

   .. code-block:: bash  

      Possible arguments are:
      -v, --version            : Software version
      -p, --port <port>        : TCP communication port with client. 
      -g, --nomodule           : [Mythen3][Gotthard2] Generic or No Module mode. 
                                 Skips detector type checks.
      -f, --phaseshift <value> : [Gotthard] only. Sets phase shift. 
      -d, --devel              : Developer mode. Skips firmware checks. 
      -u, --update             : Update mode. Skips firmware checks and initial detector setup. 
      -s, --stopserver         : Stop server. Do not use as it is created by control server 


Basics
------------

Detector Servers include:
   * Control server [default port: 1952]
      * Almost all client communication.
   * Stop server [default port: 1953]
      *  Client requests for detector status, stop acquisition, temperature, advanced read/write registers.

When using a blocking acquire command (sls_detector_acquire or Detector::acquire), the control server is blocked until end of acquisition. However, stop server commands could be used in parallel.


.. _Automatic start servers:
Automatic start 
------------------

One can start the on-board detector server automatically upon powering on the board.

#. Create a soft link to the binary on board 
   :
      .. code-block:: bash
      
         ln -sf someDetectorServervx.x.x someDetectorServer



#. Do the following depending on the detector type :

   Eiger
      .. code-block:: bash
         
         # create script in rc5.d on the board
         vi /etc/rc5.d/S50board_com.sh

         # enter the following (edit server name)
         #! /bin/sh
         /home/root/executables/eigerDetectorServer &> /dev/null &
         exit 0

   Jungfrau | Moench | CTB | Gotthard I
      .. code-block:: bash

         # Edit inittab on board
         vi /etc/inittab

         # enter the following line
         ttyS0::respawn:/./xxxDetectorServer


   Gotthard II | Mythen III
      .. code-block:: bash
         
         # create script in init.d on board
         vi /etc/init.d/S99detServer.sh

         # enter the following (edit server name)
         #! /bin/sh
         cd /root >> /dev/null
         /root/xxxDetectorServer >> /dev/null &


#. Sync, reboot and verify
   :
      .. code-block:: bash
      
         sync
         reboot

         # verify
         ps -ef | grep xxxDetectorServer
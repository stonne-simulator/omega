Hola Antonio,

aqui te resumo las ejecuciones que he realizado para la generacion de los JSON

Execution_4_VNs_32_mswitches.txt:

   -  Comando usado:  ./stonne -R=2 -S=2 -C=2 -K=32 -G=1 -N=1 -X=8 -Y=8 -T_R=2 -T_S=2 -T_C=2 -T_K=4 -T_G=1 -T_N=1 -T_X_=1 -T_Y_=1
                      -num_ms=32 -dn_bw=8 -rn_bw=8. 

   - Numero de multiplicadores: 32   

   - Numero de neuronas virtuales: 4
   
   - Notas extras: Se usan todos los multiplicadores, repartidos en 4 neuronas virutales de 8 de tamanio cada uno. 


Execution_2_VNs_32_mswitches.txt:
  
    - Comando usado: ./stonne -R=2 -S=2 -C=2 -K=32 -G=1 -N=1 -X=8 -Y=8 -T_R=2 -T_S=2 -T_C=2 -T_K=2 -T_G=1 -T_N=1 
                     -T_X_=1 -T_Y_=1 -num_ms=32 -dn_bw=8 -rn_bw=8

    - Numero de multiplicadores: 32
  
    - Numero de neuronas virtuales: 2

    - Notas extra: Hay 32 multiplicadores pero solo 16 estan en uso por medio de 2 neuronas virtuales de 8 cada una. 

Execution_2_larger_VNs_32_mswitches.txt:

     - Comando usado: /stonne -R=2 -S=2 -C=4 -K=32 -G=1 -N=1 -X=8 -Y=8 -T_R=2 -T_S=2 -T_C=4 -T_K=2 -T_G=1 
                      -T_N=1 -T_X_=1 -T_Y_=1 -num_ms=32 -dn_bw=8 -rn_bw=8

     - Numero de multiplicadores: 32

     - Numero de neuronas virtuales: 2

     - Notas extra: Esta vez solo hay 2 neuronas virtuales, pero cada una tiene tamanio 16, por lo que se usan todos los multiplicadores

Execution_2_VNs_64_mswitches.txt:

     - Comando usado: ./stonne -R=2 -S=2 -C=8 -K=32 -G=1 -N=1 -X=8 -Y=8 -T_R=2 -T_S=2 
     -T_C=8 -T_K=2 -T_G=1 -T_N=1 -T_X_=1 -T_Y_=1 -num_ms=64 -dn_bw=8 -rn_bw=8

     - Numero de multiplicadores: 64

     - Numero de neuronas virtuales: 2

     - Notas extra: Usamos 2 neuronas virtuales de 32 cada una, usando todos los multiplicadores


Execution_8_VNs_64_mswitches.txt:

     - Comando usado: ./stonne -R=2 -S=2 -C=2 -K=32 -G=1 -N=1 -X=8 -Y=8 -T_R=2 -T_S=2 
      -T_C=2 -T_K=8 -T_G=1 -T_N=1 -T_X_=1 -T_Y_=1 -num_ms=64 -dn_bw=8 -rn_bw=8

     - Numero de multiplicadores: 64

     - Numero de neuronas virtuales: 8
 
     - Notas extra: Usamos 8 neuronas virtuales de 8 multiplicadores cada una, usando asi todos los multiplicadores disponibles. 
     
Execution_2_VNs_128_mswitches.txt:

     - Comando usado: ./stonne -R=2 -S=2 -C=16 -K=32 -G=1 -N=1 -X=8 -Y=8 -T_R=2 -T_S=2 -T_C=16 
        -T_K=2 -T_G=1 -T_N=1 -T_X_=1 -T_Y_=1 -num_ms=128 -dn_bw=8 -rn_bw=8

     - Numero de multiplicadores: 128

     - Numero de neuronas virtuales: 2

     - Notas extra: Se usan dos neuronas virtuales de 64 multiplicadores cada una, usando tods los multiplicadores disponibles. 


Execution Execution_4_VNs_SIZE_27_128_mswitches.txt:

    - Comando usado: ./stonne -R=3 -S=3 -C=3 -K=32 -G=1 -N=1 -X=8 -Y=8 -T_R=3 -T_S=3 -T_C=3 
     -T_K=4 -T_G=1 -T_N=1 -T_X_=1 -T_Y_=1 -num_ms=128 -dn_bw=8 -rn_bw=8

    - Numero de multiplicadores: 128

    - Numero de neuronas virtuales: 4

    -  Notas extra: WARNING, se usan 4 neuronas virtuales de 27 multiplicadores cada una, usando 108 totales y dejando libres 20. CUIDADO, porque esta ejecucion usa los enlaces forwarding links entre los Adders de un mismo nivel. Los enlaces ESPECIALES. 




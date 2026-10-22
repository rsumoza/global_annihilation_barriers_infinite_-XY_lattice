AUDITORIA INDEPENDIENTE COMPLETA DE LA GEOMETRIA DEL SEPARADOR
21 de octubre de 2026

Objeto y alcance
---------------
check_geometry_exact.py no importa ni ejecuta el motor C++ del suplemento.
Lee los bytes binary64 little-endian como racionales diadicos exactos y usa
enteros de longitud arbitraria y Fraction de la biblioteca estandar Python.
Comprueba todas las inclusiones geometricas publicadas para los datos fijos:

  * 4 entradas de Gram y las 4 de su inversa;
  * 218 coordenadas proyectadas y 109 normas transversales al cuadrado;
  * las 380 distancias entre vecinos, sin muestreo;
  * las tres cotas escalares maximas correspondientes.

Son 715 inclusiones escalares de intervalos mas tres maximos. Se verifican
las identidades y la completitud del registro de centros/vecinos, las huellas
de todos los datos empleados y la huella congelada del manifiesto cientifico.

Identidades exactas empleadas
----------------------------
Si z es una propuesta de soporte finito extendida por cero, se calcula

  (L z)(v) = 4 z(v) - sum_{w vecino de v} z(w),
  <Dz,Dw> = sum_v z(v) (Lw)(v).

Se incluyen los cuatro grados de cada vertice, tambien en el borde del array.
En estos datos todas las propuestas se anulan en el borde; el comprobador lo
verifica explicitamente, por lo que las aristas que salen del array tienen
diferencia cero. Esta formulacion verifica las condiciones de extension y
es independiente de la acumulacion de diferencias de enlaces del motor C++.

Con q_0,q_1 las propuestas de modos, G_{ij}=<Dq_i,Dq_j>,
g(z)_i=<Dz,Dq_i>, a(z)=G^{-1}g(z), se calcula

  ||Yz||^2 = ||Dz||^2 - g(z)^T G^{-1} g(z).

Para los vecinos se calcula primero exactamente

  ||D(z_i-z_j)||^2 = ||Dz_i||^2 + ||Dz_j||^2 - 2<Dz_i,Dz_j>,

y se resta la componente proyectada. No se aproxima ninguna raiz cuadrada:
una cota decimal [l,h] para una distancia se acepta solo si
0 <= l <= h y l^2 <= distancia^2 <= h^2, en aritmetica racional exacta.

Ejecucion
---------
Desde este directorio, para la salida historica incluida:

python3 check_geometry_exact.py --package ../xy_global_separator \
  --geometry ../xy_global_separator/results_provided/geometry.json \
  --output exact_geometry_archived.json

Para la salida reconstruida por el supervisor en esta auditoria:

python3 check_geometry_exact.py --package ../xy_global_separator \
  --geometry ../fresh_results/geometry.json \
  --output exact_geometry_fresh.json

La ejecucion historica paso todas las inclusiones en aproximadamente 3.4 s.
La ejecucion sobre geometry.json reconstruido el 21 de octubre de 2026
tambien paso las 715 inclusiones y los tres maximos (aproximadamente 2.9 s).
Las salidas geometricas historica y fresca coinciden byte por byte:
SHA256 13076bc4aaeebcdda84c5e4aa795a47e0cfd5fc992e7acaa72414f8f85f9a186.
El archivo JSON identifica cada dato y el codigo del comprobador mediante SHA256.
Los controles usan excepciones explicitas; Python -O no los elimina.

Controles adversariales
----------------------
check_negative_controls.py modifica en archivos temporales una coordenada
para sustituir su intervalo por [0,0], que excluye su valor exacto, y elimina
por separado un vecino del registro. Ambos casos se ejecutan con Python -O:
deben fallar por la inclusion incorrecta y por el registro incompleto,
respectivamente, sin producir un informe de aprobacion. La ejecucion dio
el rechazo esperado para ambos controles, tanto con la salida historica
como con la salida fresca. Para reproducir:

python3 check_negative_controls.py --package ../xy_global_separator \
  --geometry ../xy_global_separator/results_provided/geometry.json \
  --output negative_controls.json

Limitacion esencial
------------------
Esto es una verificacion independiente completa de la geometria finita
certificada por geometry.json, incluyendo su significado en la norma de
Dirichlet de la red infinita por extension cero. No verifica por si sola
las cotas analiticas del minimo infinito, las cartas energeticas, el cierre
topologico del separador ni el camino superior de aniquilacion. No constituye
una formalizacion en asistente de pruebas.

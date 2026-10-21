# Setup en calafate (Linux)

Procedimiento para dejar el proyecto compilando y resolviendo FMO en
`calafate.inf.ucv.cl`. Todo lo que sigue fue ejecutado y verificado; los valores
concretos (versiones, rutas, tamaños) corresponden a esa ejecución.

El repositorio vive en `~/emili_imrt/EMILI-IMRT`.

## 1. Por qué no alcanza con clonar

Hay cuatro diferencias entre el entorno de desarrollo (macOS ARM) y calafate:

| # | Diferencia | Impacto |
|---|---|---|
| 1 | `libampl.dylib` es Mach-O | Linux necesita el ELF `libampl.so` |
| 2 | `imrt_fmo.cpp` asume `python3.9` en el fallback | calafate tiene Python 3.6.8 |
| 3 | `ampl_gurobi/.venv/` está en `.gitignore` | No viaja con el repo |
| 4 | `ampl.lic` está en `.gitignore` | No viaja con el repo |

La diferencia 1 se resolvió en el repositorio: `CMakeLists.txt` ahora elige la
librería por plataforma y ambas están vendorizadas bajo
`third_party/ampl_cppapi/lib/`. Las diferencias 2, 3 y 4 se resuelven con este
procedimiento.

## 2. Instancias

`instances/` está en `.gitignore` por tamaño (11 GB, 2174 archivos). Se copia
por `rsync` desde la máquina de desarrollo:

```bash
rsync -avhP -z \
  instances/CERR_Prostate/ \
  calafate:emili_imrt/EMILI-IMRT/instances/CERR_Prostate/
```

`-z` no es opcional en la práctica: el contenido es texto plano
(`voxel  beamlet  dosis`) y comprime 5.4x, lo que bajó la transferencia de
~18 minutos a 3:20.

## 3. Módulos AMPL

Los módulos no están en PyPI estándar sino en el índice propio de AMPL. Se
descargan **desde la máquina de desarrollo** y se transfieren, de modo que en
calafate no se descarga nada:

```bash
pip3 download --index-url https://pypi.ampl.com \
  --extra-index-url https://pypi.org/simple \
  --only-binary=:all: --platform manylinux2014_x86_64 \
  --python-version 3.6 --implementation py --no-deps \
  ampl_module_base ampl_module_gurobi -d ./ampl_linux

rsync -avh -z ./ampl_linux/*.whl calafate:ampl_setup/
```

Total: ~18 MB. Los wheels tienen tag `py3`, por lo que Python 3.6 los acepta.

`libampl.so` sale del wheel de `amplpy` para Linux, en
`amplpy/amplpython/cppinterface/lib/amd64/libampl.so`. Ya está vendorizado en el
repositorio, así que este paso solo es necesario al actualizar la versión de
`amplpy`. Si se actualiza, deben actualizarse `.so`, `.dylib` y los headers
juntos: provienen del mismo wheel y comparten ABI.

## 4. Entorno virtual

```bash
cd ~/emili_imrt/EMILI-IMRT
python3 -m venv ampl_gurobi/.venv
ampl_gurobi/.venv/bin/pip install --no-index \
  ~/ampl_setup/ampl_module_base-*.whl \
  ~/ampl_setup/ampl_module_gurobi-*.whl
```

Queda instalado bajo `ampl_gurobi/.venv/lib/python3.6/site-packages/`.

## 5. Licencia

La licencia es el bundle académico del curso Optimización 2 ICD
(Bundle #7591.8177, vence 2027-01-16). Es un archivo local: no consulta
servidor de tokens ni tiene límite de usuarios concurrentes.

```bash
cp ampl.lic ampl_gurobi/.venv/lib/python3.6/site-packages/ampl_module_base/bin/
chmod 600 ampl_gurobi/.venv/lib/python3.6/site-packages/ampl_module_base/bin/ampl.lic
```

El `chmod 600` es deliberado: calafate es multiusuario.

**No usar `/opt/gurobi1001`.** Esa instalación (Gurobi 10.0.1) pertenece a un
uid sin entrada en `passwd`, no tiene licencia y falla con Error 10009. El
Gurobi que usa este proyecto es el 13.x que trae `ampl_module_gurobi` dentro del
venv, habilitado por el bundle de AMPL. No se debe exportar `GRB_LICENSE_FILE`
ni agregar `/opt/gurobi1001/linux64/lib` al `LD_LIBRARY_PATH`.

## 6. Compilación

CMake en calafate es 3.11.4, que no soporta los flags `-S` / `-B`
(introducidos en 3.13). Hay que usar la forma clásica:

```bash
cd ~/emili_imrt/EMILI-IMRT
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```

Verificar que el binario resolvió la librería correcta:

```bash
ldd build/emili | grep ampl
# libampl.so => .../third_party/ampl_cppapi/lib/libampl.so
```

## 7. Variables de entorno

`imrt_fmo.cpp` tiene como fallback una ruta con `python3.9`, que en calafate no
existe. Las variables `EMILI_AMPL_BIN_DIR` y `EMILI_GUROBI_BIN` son la costura
prevista para eso y evitan tener que tocar el código:

```bash
V=~/emili_imrt/EMILI-IMRT/ampl_gurobi/.venv/lib/python3.6/site-packages
export EMILI_AMPL_BIN_DIR=$V/ampl_module_base/bin
export EMILI_GUROBI_BIN=$V/ampl_module_gurobi/bin/gurobi
```

Conviene agregarlas al `~/.bashrc` o al script de lanzamiento de experimentos.
Sin ellas el binario compila y corre, pero falla al inicializar AMPL.

## 8. Verificación

```bash
cd ~/emili_imrt/EMILI-IMRT
./build/emili instances/CERR_Prostate baoimrt 4 csv /tmp/smoke.csv \
  ils first irandomk locmin nangshift 5 tmaxiter 1 \
  prangshift 5 3 baoimprove rejectrepeated rnds 1
```

Debe terminar con código 0 y emitir la tabla de restricciones DVH y los índices
de calidad. Con `tmaxiter 1` el plan resultante es malo y las restricciones
aparecen como `VIOL`: eso es esperado. Lo que se está verificando es que los FMO
se resuelven y que el reporte se genera, no la calidad del plan.



# Curvatura de paraboloides

Este es un fork de [Pececillo](https://gitlab.com/truchas/pececillo) que reproduce los resultados del artículo "A Paraboloid Fitting Technique for Calculating Curvature from Piecewise-Linear Interface Reconstructions on 3D Unstructured Meshes" de Z. Jibben, N. N. Carlson y M. M. Francois, 2017.

## Compilación

Depende de las bibliotecas necesarias para compilar Pececillo. Puede encontrar una superbuild [aquí](https://gitlab.com/truchas/pececillo-tpl/tags/v1). Las pruebas se ejecutan correctamente con los compiladores de Intel v17.0.4.

```bash
$ mkdir build
$ cd build
$ cmake -DCMAKE_PREFIX_PATH=/your/pececillo-tpl/install/path ..
$ make
```

## Ejecución

```bash
$ cd run
$ ../build/src/parabfit/pececillo-parabfit
```

De forma predeterminada, esto calculará las normas de error de curvatura para el caso de prueba del elipsoide en la malla tetraédrica más gruesa, lo que debería tardar entre 5 y 20 minutos, dependiendo de su hardware.

Para reproducir todos los resultados del artículo, primero debe descargar todas las mallas tetraédricas y hexaédricas deformadas [aquí](https://zenodo.org/record/1135295/files/meshes.tar.gz) y descomprimirlas en el directorio `run`. A continuación, ejecute `$ ../build/src/parabfit/pececillo-parabfit full`. Esto tardará entre días y semanas en calcularse en un solo puesto de trabajo, debido al alto costo de inicializar los campos de fracción de volumen y, en menor medida, de calcular los vectores normales mediante LVIRA. Puede acelerar este proceso descargando en su lugar los campos de fracción de volumen y los vectores normales [aquí](https://zenodo.org/record/1135295/files/fields.tar.gz) y descomprimiéndolos en el directorio `run`. Al cargar estos campos, todos los resultados del artículo se calcularán en unas pocas horas.

Tanto los archivos `tarball` de campos como de mallas están identificados por [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1135295.svg)](https://doi.org/10.5281/zenodo.1135295).

## Algoritmo de ajuste de paraboloides en el código fuente

Para encontrar el código fuente del algoritmo de ajuste de paraboloides, `consulte lo siguiente en el directorio `src/lib`:`

- `surface_type.F90`, función `local_patch` devuelve una colección de polígonos de reconstrucción de interfaz que rodean una celda objetivo.
- `interface_patch_type.F90`, función `curvature_from_patch` toma una colección de polígonos y `gestiona otras rutinas para ajustar un paraboloides y calcular la curvatura.
- `paraboloid_type.F90`, subrutina `volumetricFit` ajusta un paraboloides a un conjunto dado de polígonos.
- `paraboloid_type.F90`, función `curvature` calcula la curvatura en un punto del paraboloides.

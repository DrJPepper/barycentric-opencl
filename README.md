# Barycentric-OpenCL

README forthcoming

## Features

## Installation

## Usage

`cmake -DDISABLE_DEMO=ON ..`

## TODO

* ☐ Improve memory management and resource cleanup
* ☐ Check available VRAM and determine max settings
* ☐ Allow user to request either max memory usage or max concurrent points
* ☐ Run optimization check to see at what point it becomes worth running on OpenCL vs CPU
* ☐ Auto switch between CPU, OCL and multiple iterations of OCL
* ☐ Decouple barycentric code further by allowing user to pass in their own vertices and faces instead of only working with input models
* ☐ Stop just returning 0,0/0,0,0 when barycentric calculations fail/error
* ☐ Fall back to CPU when a point fails on OCL and print a more detailed error message
<!--☑-->

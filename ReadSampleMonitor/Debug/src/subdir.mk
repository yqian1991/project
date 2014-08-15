################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/dataPreProcess.c \
../src/genePreProcess.c \
../src/monitor.c \
../src/test.c 

OBJS += \
./src/dataPreProcess.o \
./src/genePreProcess.o \
./src/monitor.o \
./src/test.o 

C_DEPS += \
./src/dataPreProcess.d \
./src/genePreProcess.d \
./src/monitor.d \
./src/test.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



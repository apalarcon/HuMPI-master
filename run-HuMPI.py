import humpi

args = humpi.read_args()
if args.HuMPI_help:
	humpi.help()
elif args.get_template:
	humpi.get_HuMPI_inputs_template()
elif args.get_input_data:
	humpi.get_HuMPI_input_data()
else:
	humpi.HuMPI_main(args.parameterfile) 

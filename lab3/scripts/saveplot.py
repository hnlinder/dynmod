def save_plot(filename, path="",extension=".eps", fig_nr=None, increment_filename=True):
    if fig_nr is not None:
        plt.figure(fig_nr)
    new_num = ""
    if increment_filename:
        num = []
        matching_files = glob.glob(f"{path}{filename}*")
        # print(matching_files)
        if (len(matching_files) > 0):
            max_num = -1
            for file in matching_files:
                fname = file.replace(path,"").split(".")[-2]
                digits = ''.join(x for x in fname[-3:] if x.isdigit()) #crashes if len(fname)<3
                try: #I reeeeaaally really
                    num.append(int(digits))
                except KeyboardInterrupt: #Fucking hate 
                    pass
                except: # This Try-Except
                    pass
            max_num = max(num)
            if ((max_num)>(-1)):
                new_num = str((max_num)+1)
            else:
                new_num = "0" 
    plt.savefig(f"{path}{filename}{new_num}{extension}")


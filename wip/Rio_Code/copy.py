axs[i,j].set_xlabel("Time [h]",   fontdict={'family':'DejaVu Sans','size':fs})
        axs[i,j].set_ylabel("Amount [mol]",   fontdict={'family':'DejaVu Sans','size':fs})
        labels = axs[i,j].get_xticklabels() + axs[i,j].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
        axs[i,j].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
from fengalotl.fct.load import get_data
import matplotlib.pyplot as plt
import scanpy as sc

def plot_gene_expression(input):  

    adata = get_data(input)
    plot_genes = list(input.select_gene_expression())

    if adata is None or not plot_genes:
        return

    plot_ex = sc.pl.dotplot(adata,
                            var_names = plot_genes,
                            swap_axes = True,
                            groupby = input.select_resolution(),
                            return_fig = True
                            )

    main_ax_ex = plot_ex.get_axes()['mainplot_ax']
    y_axis_labels = [tick.get_text() for tick in main_ax_ex.get_yticklabels()]
    main_ax_ex.set_yticklabels(y_axis_labels)

    plt.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.95)

def plot_de(input):  

    if not (adata := get_data(input)):
        return None
    
    key = 'rank_' + input.select_resolution()
    if key not in adata.uns:
        # Run differential expression if not already computed
        sc.tl.rank_genes_groups(adata, groupby=input.select_resolution(), key_added=key)
    
    plot_de = sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes = input.slider_n_genes(),
        key = key,
        min_logfoldchange=input.slider_lfc(),
        return_fig = True)

    main_ax_de = plot_de.get_axes()['mainplot_ax']
    x_axis_labels = [tick.get_text() for tick in main_ax_de.get_xticklabels()]
    main_ax_de.set_xticklabels(x_axis_labels)
    
    plt.subplots_adjust(top=1, bottom=0.2, left=0.05, right=0.95)

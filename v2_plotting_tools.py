import matplotlib.pyplot as plt
import numpy as np
import math

def heat_map(data, colormap, ** keywords):
    # keywords
        # colormap
        # Title
        # X_title
        # Y_title
        # X_labels
        # Y_labels
        # FigSize
        return None

def scatter_plot(x_data, y_data, **keywords): #expects list for y_data
    '''keywords:
        Colors = 'rainbow' 'ACGT' 'AGTN'
        Title
        X_title
        Y_title
        FigSize
        Y_log
    '''
    # colors
    Colors = keywords.pop('Colors','rainbow')
    if Colors == 'rainbow':
        color_list = ['xkcd:red','xkcd:orange','xkcd:yellow',
                      'xkcd:green','xkcd:blue','xkcd:purple']
        Colors = color_list[:len(y_data)]
    elif Colors == 'ACGT':
        Colors = ['g','b','0','r']
    elif Colors == 'ACGTN':
        Colors = ['g','b','0','r','0.6']


    # plotting keywords
    FigSize = keywords.pop('FigSize',None)
    Marker = keywords.pop('Marker','o')
    MarkerSize = keywords.pop('MarkerSize',5)
    LineStyle = keywords.pop('LineStyle','')
    Y_log = keywords.pop('Y_log',False)

    # plot each series, find bounds
    y_min, y_max = 1, 0
    x_min, x_max = min(x_data), max(x_data)

    fig, ax = plt.subplots(figsize = FigSize)
    for y_series, color in zip (y_data, Colors):
        y_series_min, y_series_max = min(y_series), max(y_series)
        if y_series_min < y_min:
            y_min = y_series_min
        if y_series_max > y_max:
            y_max = y_series_max
        ax.plot(x_data, y_series, marker=Marker, markersize=MarkerSize, markerfacecolor=color,
                linestyle=LineStyle, markeredgewidth=0)

    if Y_log:
        lower_bound_y = 10**(int(math.log10(y_min))-1)
        upper_bound_y = 10**(int(math.log10(y_max)))
        ax.set_yscale('log')
    else:
        lower_bound_y = int(y_min) - int(y_min)*0.05
        upper_bound_y = int(y_max)*1.05
        if upper_bound_y < y_max:
            upper_bound_y = int(y_max)+1
        if lower_bound_y == 0:
            lower_bound_y = upper_bound_y*-0.05

    lower_bound_x = int(x_min) - int(x_min)*0.05
    upper_bound_x = int(x_max) + int(x_max)*0.05
    if upper_bound_x < x_max:
        upper_bound_x = int(x_max)+1
    if lower_bound_x == 0:
        lower_bound_x = upper_bound_x*-0.05

    ax.set_ylim(lower_bound_y, upper_bound_y)
    ax.set_xlim(lower_bound_x, upper_bound_x)


    # axis labels and titles
    Title = keywords.pop('Title',None)
    if Title:
        ax.set_title(Title)
    X_title = keywords.pop('X_title',None)
    if X_title:
        ax.set_xlabel(X_title)
    Y_title = keywords.pop('Y_title',None)
    if Y_title:
        ax.set_ylabel(Y_title)

    plt.tight_layout()
    return fig, ax


def bar_chart(category_labels, data_series, ** keywords):
    '''keywords:
        Colors = 'rainbow' 'ACGT' 'AGTN'
        Title
        X_title
        Y_title
        FigSize
        Y_log
    '''
    # colors
    Colors = keywords.pop('Colors','rainbow')
    if Colors == 'rainbow':
        color_list = ['xkcd:red','xkcd:orange','xkcd:yellow',
                      'xkcd:green','xkcd:blue','xkcd:purple']
        Colors = color_list[:len(data_series)]
    elif Colors == 'ACGT':
        Colors = ['g','b','0.2','r']
    elif Colors == 'ACGTN':
        Colors = ['g','b','0.2','r','0.6']

    # plotting keywords
    FigSize = keywords.pop('FigSize',None)
    Marker = keywords.pop('Marker','o')
    MarkerSize = keywords.pop('MarkerSize',5)
    LineStyle = keywords.pop('LineStyle','')
    Y_log = keywords.pop('Y_log',False)

    # plot each series
    y_max = 0
    x_pos = np.arange(len(data_series[0]))
    series_count = len(data_series) + 1
    x_pos_subs = np.arange(0,1,1/series_count)

    fig, ax = plt.subplots()
    for series, Color, x_pos_sub in zip(data_series, Colors, x_pos_subs):
        plt.bar(x_pos + x_pos_sub, series, color=Color, edgecolor='0',
                width=(1/len(x_pos_subs)), align='edge')
        if max(series) > y_max:
            y_max = max(series)

    if Y_log:
        ax.set_yscale('log')

    # axis labels and titles
    Title = keywords.pop('Title',None)
    if Title:
        ax.set_title(Title)
    X_title = keywords.pop('X_title',None)
    if X_title:
        ax.set_xlabel(X_title)
    Y_title = keywords.pop('Y_title',None)
    if Y_title:
        ax.set_ylabel(Y_title)

    # x axis labels
    x_label_pos = []
    for x in x_pos:
        x_label_pos.append(x+0.5-(1/len(x_pos_subs)/2))
    ax.set_xticks(x_label_pos)
    ax.set_xticklabels(category_labels)

    # add error bar capability, log capability

    plt.tight_layout()

    return fig, ax


def stacked_bar_chart(category_labels, data_series, ** keywords):
    '''keywords:
        Colors = 'rainbow' 'ACGT' 'AGTN'
        Title
        X_title
        Y_title
        FigSize
        Y_log
    '''
    # colors
    Colors = keywords.pop('Colors','rainbow')
    if Colors == 'rainbow':
        color_list = ['xkcd:red','xkcd:orange','xkcd:yellow',
                      'xkcd:green','xkcd:blue','xkcd:purple']
        Colors = color_list[:len(data_series)]
    elif Colors == 'ACGT':
        Colors = ['g','b','0.2','r']
    elif Colors == 'ACGTN':
        Colors = ['g','b','0.2','r','0.6']

    # plotting keywords
    FigSize = keywords.pop('FigSize',None)
    Marker = keywords.pop('Marker','o')
    MarkerSize = keywords.pop('MarkerSize',5)
    LineStyle = keywords.pop('LineStyle','')
    Y_log = keywords.pop('Y_log',False)

    # plot each series
    x_pos = np.arange(len(data_series[0]))
    bottoms = list(np.zeros(len(data_series[0])))

    fig, ax = plt.subplots()
    for series, Color in zip(data_series, Colors):
        plt.bar(x_pos, series, color=Color, edgecolor='0', bottom=bottoms)
        for d in range(len(series)):
            bottoms[d] += series[d]

    # axis labels and titles
    Title = keywords.pop('Title',None)
    if Title:
        ax.set_title(Title)
    X_title = keywords.pop('X_title',None)
    if X_title:
        ax.set_xlabel(X_title)
    Y_title = keywords.pop('Y_title',None)
    if Y_title:
        ax.set_ylabel(Y_title)

    # x axis labels
    ax.set_xticks(x_pos)
    ax.set_xticklabels(category_labels)

    plt.tight_layout()

    return fig, ax

if __name__ == "__main__":
    test_x_data = [10,20,30]
    test_y_data = [[1,2,3],[4,5,6],[7,7,7],[1,1,6]]


    fig, ax = scatter_plot(test_x_data, test_y_data, FigSize=(4, 2), Y_log=False)


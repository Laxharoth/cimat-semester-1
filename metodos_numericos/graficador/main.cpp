#include "point.hpp"
#include "simple_function_analizer.hpp"
#include "canvas_facade.hpp"

#include <array>
#include <iostream>
#include <gtk/gtkx.h>
#include <gtk/gtk.h>

// TODO: ADD scale to function to fill heigh and width
// TODO: ADD x and y axis
// TODO: ADD marks in x and y axis showing the value

using canvas_facade::CanvasFacade;

void init_points(GtkWidget *graph);
gboolean on_graph_draw(GtkWidget *graph, cairo_t *cr, gpointer data);
void on_destroy();
CanvasFacade drawer;

std::array<point, PARTITIONS_NUM> points;


int main(int argc, char *argv[]){
    gtk_init(&argc, &argv);

    GtkBuilder *builder = gtk_builder_new_from_file("graph.glade");
    GtkWidget  *window = GTK_WIDGET(gtk_builder_get_object(builder, "window"));
    GtkWidget  *graph = GTK_WIDGET(gtk_builder_get_object(builder, "graph"));

    g_object_unref(builder);

    g_signal_connect(window, "destroy", G_CALLBACK(on_destroy), NULL);
    g_signal_connect(graph, "draw", G_CALLBACK(on_graph_draw), NULL);

    gtk_window_set_keep_above( GTK_WINDOW(window), TRUE );

    gtk_widget_show(window);
    init_points(graph);
    gtk_widget_queue_draw(graph);

    gtk_main();

    return EXIT_SUCCESS;
}

void init_points(GtkWidget *graph){
    double min_y{};
    double max_y{};
    gint width{}, height{};
    width = gtk_widget_get_allocated_width(graph);
    height= gtk_widget_get_allocated_height(graph);
    analize_single_var_function([](double x){ return x*x; }, -10, 10, min_y, max_y, points);
    drawer.min_x = points[0].x;
    drawer.min_y = min_y;
    drawer.set_origin( canvas_facade::choose_origin_helper(points[0].x, min_y, height) );
}

void on_destroy(){
    gtk_main_quit();
}
gboolean on_graph_draw(GtkWidget *graph, cairo_t *cr, gpointer data){
    g_printerr("origin: (%f,%f)" ,drawer.origin.x, drawer.origin.y);
    for(size_t i = 0; i < points.size() - 1 ; ++i){
        drawer.draw_line(cr, points[i], points[i +1 ]);
    }
    return FALSE;
}

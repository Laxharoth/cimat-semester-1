#include "point.hpp"
#include "simple_function_analizer.hpp"
#include "canvas_facade.hpp"

#include <array>
#include <iostream>
#include <gtk/gtkx.h>
#include <gtk/gtk.h>

using canvas_facade::CanvasFacade;

void init_points(GtkWidget *graph, double (*single_variable_func)(double), const double &min_x, const double &max_x);
gboolean on_graph_draw(GtkWidget *graph, cairo_t *cr, gpointer data);
std::string trim_double_to_str(const double &num, const int &precision);
void on_destroy();

CanvasFacade drawer;
std::array<point, PARTITIONS_NUM> points;

gint main(int argc, char *argv[]){
    gtk_init(&argc, &argv);
    GtkBuilder *builder = gtk_builder_new_from_file("graph.glade");
    GtkWidget  *window = GTK_WIDGET(gtk_builder_get_object(builder, "window"));
    GtkWidget  *graph = GTK_WIDGET(gtk_builder_get_object(builder, "graph"));
    g_object_unref(builder);
    g_signal_connect(window, "destroy", G_CALLBACK(on_destroy), NULL);
    g_signal_connect(graph, "draw", G_CALLBACK(on_graph_draw), NULL);
    gtk_window_set_keep_above( GTK_WINDOW(window), TRUE );
    gtk_widget_show(window);

    init_points(graph,[](double x){ return 1/x; }, -1.0, 1.0);

    gtk_main();

    return EXIT_SUCCESS;
}
void draw_axis(cairo_t * cr, const gint width, const gint height, const double &min_x, const double &max_x, const double &min_y, const double max_y){
    const gint axis_lines = 10;
    const gint middle_width   = width / 2, middle_height = height/ 2;
    const gint distance_width = width / axis_lines, distance_height = height/ axis_lines;
    const double middle_x   = (max_x + min_x) / 2, middle_y = (max_y + min_y)/ 2;
    const double change_x   = (max_x - min_x) / axis_lines, change_y = (max_y - min_y) / axis_lines;
    cairo_move_to(cr, 0, middle_height); cairo_line_to(cr, width , middle_height);
    cairo_move_to(cr, middle_width, 0); cairo_line_to(cr, middle_width , height);
    for(gint line = 0; line < axis_lines / 2 ; ++line ){
        // right
        cairo_move_to(cr, middle_width + distance_width * line, middle_height - 4);
        cairo_line_to(cr, middle_width + distance_width * line, middle_height + 4);
        cairo_move_to(cr, middle_width + distance_width * line + 4, middle_height + 10);
        cairo_show_text(cr, trim_double_to_str(middle_x + change_x * line, 2).c_str());
        // up
        cairo_move_to(cr, middle_width - 4, middle_height - distance_height * line);
        cairo_line_to(cr, middle_width + 4, middle_height - distance_height * line);
        cairo_move_to(cr, middle_width + 4, middle_height - distance_height * line - 10);
        cairo_show_text(cr, trim_double_to_str(middle_y + change_y * line, 2).c_str());
        if(line == 0) continue;
        // left
        cairo_move_to(cr, middle_width - distance_width * line, middle_height - 4);
        cairo_line_to(cr, middle_width - distance_width * line, middle_height + 4);
        cairo_move_to(cr, middle_width - distance_width * line + 4, middle_height + 10);
        cairo_show_text(cr, trim_double_to_str(middle_x - change_x * line, 2).c_str());
        // down
        cairo_move_to(cr, middle_width - 4, middle_height + distance_height * line);
        cairo_line_to(cr, middle_width + 4, middle_height + distance_height * line);
        cairo_move_to(cr, middle_width + 4, middle_height + distance_height * line - 10);
        cairo_show_text(cr, trim_double_to_str(middle_y - change_y * line, 2).c_str());
    }
    cairo_stroke(cr);
}

void init_points(GtkWidget *graph, double (*single_variable_func)(double), const double &min_x, const double &max_x){
    double min_y{}, max_y{};
    analize_single_var_function(single_variable_func, min_x, max_x, min_y, max_y, points);
    drawer.min_value.x = min_x; drawer.min_value.y = min_y;
    drawer.max_value.x = max_x; drawer.max_value.y = max_y;
    gtk_widget_queue_draw(graph);
}

void on_destroy(){
    gtk_main_quit();
}

gboolean on_graph_draw(GtkWidget *graph, cairo_t *cr, gpointer data){
    gint width = gtk_widget_get_allocated_width(graph);
    gint height= gtk_widget_get_allocated_height(graph);
    drawer.scale_x = canvas_facade::choose_scale_helper(drawer.min_value.x, drawer.max_value.x, width);
    drawer.scale_y = canvas_facade::choose_scale_helper(drawer.min_value.y, drawer.max_value.y, height);
    drawer.set_origin( canvas_facade::choose_origin_helper(drawer.min_value.x, drawer.min_value.y,drawer.scale_x,drawer.scale_y, height) );
    drawer.draw_multiline((void*)cr, points.data(), points.size());
    draw_axis(cr, width, height, drawer.min_value.x, drawer.max_value.x, drawer.min_value.y, drawer.max_value.y);
    return FALSE;
}

std::string trim_double_to_str(const double &num, const int &precision){
    std::string n_to_s = std::to_string(num);
    int cut =  n_to_s.find(".");
    if( cut < 0 || cut + precision + 1 > n_to_s.length() )return n_to_s;
    return n_to_s.substr(0, cut + precision + 1);
}

#include "point.hpp"
#include "simple_function_analizer.hpp"
#include "canvas_facade.hpp"
#include "function_wrapper/function_wrapper.hpp"
#include "function_wrapper/parser_wrapper.hpp"
#include "fparser/fparser.hh"

#include <array>
#include <iostream>
#include <gtk/gtkx.h>
#include <gtk/gtk.h>

#include <cmath>
#define Tolerancia_y 10E-6

using canvas_facade::CanvasFacade;

/**
 * @brief Mapea los puntos en el intervalo especificado y actualiza los minimos y maximos de x y f(x)
 */
void init_points(GtkWidget *graph, double (*single_variable_func)(double), double min_x, double max_x);
/**
 * @brief Dibuja en la superficie cairo
 */
gboolean on_graph_draw(GtkWidget *graph, cairo_t *cr, gpointer data);
/**
 * @brief vuelve a dibujar la funcion
 */
void on_btn_graph_click(GtkWidget *btn, gpointer data);
/**
 * @brief convierte un double a string con n puntos decimales (no redondea)
 */
std::string trim_double_to_str(const double &num, const int &precision);
/**
 * @brief *cambia la funcion que se muestra
 */
void on_radio_original_toggle(GtkWidget *widget, gpointer data);
void on_destroy();

CanvasFacade drawer;
FunctionParser parser;
FunctionWrapper *current_function;
FunctionParserAdapter_to_FWrapper original(&parser);
Derivative derivative(&original);
Derivative derivative_2(&derivative);

std::array<point, PARTITIONS_NUM+1> points;
GtkWidget  *graph;
GtkWidget  *btn_graph;
GtkWidget  *min_x_limit;
GtkWidget  *max_x_limit;
GtkWidget  *function_rep;
GtkWidget  *radio_original;
GtkWidget  *radio_derivada;
GtkWidget  *radio_derivada_2;

auto func = [](double x){ return 1/x; };

gint main(int argc, char *argv[]){
    current_function = &original;
    gtk_init(&argc, &argv);
    GtkBuilder *builder = gtk_builder_new_from_file("graph.glade");
    GtkWidget  *window = GTK_WIDGET(gtk_builder_get_object(builder, "window"));
    graph = GTK_WIDGET(gtk_builder_get_object(builder, "graph"));
    min_x_limit = GTK_WIDGET(gtk_builder_get_object(builder, "min_x_limit"));
    max_x_limit = GTK_WIDGET(gtk_builder_get_object(builder, "max_x_limit"));
    function_rep = GTK_WIDGET(gtk_builder_get_object(builder, "function_rep"));
    radio_original = GTK_WIDGET(gtk_builder_get_object(builder, "radio_original"));
    radio_derivada = GTK_WIDGET(gtk_builder_get_object(builder, "radio_derivada"));
    radio_derivada_2 = GTK_WIDGET(gtk_builder_get_object(builder, "radio_derivada_2"));
    btn_graph = GTK_WIDGET(gtk_builder_get_object(builder, "btn_graph"));
    g_object_unref(builder);
    g_signal_connect(window, "destroy", G_CALLBACK(on_destroy), NULL);
    g_signal_connect(graph, "draw", G_CALLBACK(on_graph_draw), NULL);
    g_signal_connect(btn_graph, "clicked", G_CALLBACK(on_btn_graph_click), NULL);
    g_signal_connect(radio_original, "toggled", G_CALLBACK(on_radio_original_toggle), (gpointer)"1");
    g_signal_connect(radio_derivada, "toggled", G_CALLBACK(on_radio_original_toggle), (gpointer)"2");
    g_signal_connect(radio_derivada_2, "toggled", G_CALLBACK(on_radio_original_toggle), (gpointer)"3");
    gtk_window_set_keep_above( GTK_WINDOW(window), FALSE );
    gtk_widget_show(window);

    //dibuja la grafica por primera vez
    on_btn_graph_click( btn_graph, NULL );

    gtk_main();

    return EXIT_SUCCESS;
}
void draw_axis(cairo_t * cr, const gint width, const gint height, const double &min_x, const double &max_x, const double &min_y, const double max_y){
    const gint axis_lines = 10;
    const gint middle_width   = width / 2, middle_height = height/ 2;
    const gint distance_width = width / axis_lines, distance_height = height/ axis_lines;
    const double middle_x   = (max_x + min_x) / 2, middle_y = (max_y + min_y)/ 2;
    const double change_x   = (max_x - min_x) / axis_lines, change_y = (max_y - min_y) / axis_lines;
    cairo_set_line_width(cr, 0.5);
    if(middle_x == 0 && middle_y == 0)
        cairo_set_line_width(cr, 1.5);
    cairo_set_source_rgb(cr,0.0,0.0,0.0);
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
    cairo_set_line_width(cr, 0.5);
    for(gint line = 0; line < axis_lines / 2 ; ++line ){
        if(line == 0) continue;
        // right
        cairo_move_to(cr, middle_width + distance_width * line, 0);
        cairo_line_to(cr, middle_width + distance_width * line, height);
        // up
        cairo_move_to(cr, 0, middle_height - distance_height * line);
        cairo_line_to(cr, width + 4, middle_height - distance_height * line);
        // left
        cairo_move_to(cr, middle_width - distance_width * line, 0);
        cairo_line_to(cr, middle_width - distance_width * line, height);
        // down
        cairo_move_to(cr, 0, middle_height + distance_height * line);
        cairo_line_to(cr, width, middle_height + distance_height * line);
    }
    cairo_stroke(cr);
}

void init_points(GtkWidget *graph, double (*single_variable_func)(double), double min_x, double max_x){
    try{
        parser.Parse(gtk_entry_get_text (GTK_ENTRY (function_rep)), "x");
    }catch(...){
        parser.Parse("x", "x");
    }
    if(min_x > max_x) std::swap(min_x, max_x);
    drawer.min_value.x = min_x; 
    drawer.max_value.x = max_x;
    analize_single_var_function([](double x){return current_function->eval(x);}, min_x, max_x, drawer.min_value.y, drawer.max_value.y, points);
    if(drawer.max_value.y - drawer.min_value.y <= Tolerancia_y){
        drawer.min_value.y -=1;
        drawer.max_value.y +=1;
    }
    gtk_widget_queue_draw(graph);
}

void on_destroy(){
    gtk_main_quit();
}

gboolean on_graph_draw(GtkWidget *graph, cairo_t *cr, gpointer data){
    const gint width = gtk_widget_get_allocated_width(graph);
    const gint height= gtk_widget_get_allocated_height(graph);
    drawer.scale_x = canvas_facade::choose_scale_helper(drawer.min_value.x, drawer.max_value.x, width);
    drawer.scale_y = canvas_facade::choose_scale_helper(drawer.min_value.y, drawer.max_value.y, height);
    drawer.set_origin( canvas_facade::choose_origin_helper(drawer.min_value.x, drawer.min_value.y,drawer.scale_x,drawer.scale_y, height) );
    drawer.draw_multiline((void*)cr, points.data(), points.size());
    draw_axis(cr, width, height, drawer.min_value.x, drawer.max_value.x, drawer.min_value.y, drawer.max_value.y);
    return FALSE;
}
void on_btn_graph_click(GtkWidget *btn, gpointer data){
    const double min_x = gtk_spin_button_get_value (GTK_SPIN_BUTTON (min_x_limit));
    const double max_x = gtk_spin_button_get_value (GTK_SPIN_BUTTON (max_x_limit));
    try{
        if(min_x == max_x)return;
        init_points(graph, func, min_x, max_x);
    }catch(std::exception& e){}
}
void on_radio_original_toggle(GtkWidget *widget, gpointer data){
    if(data == "1"){ current_function = &original; }
    if(data == "2"){ current_function = &derivative; }
    if(data == "3"){ current_function = &derivative_2; }
    on_btn_graph_click( btn_graph, NULL );
}
std::string trim_double_to_str(const double &num, const int &precision){
    std::string n_to_s = std::to_string(num);
    int cut =  n_to_s.find(".");
    if( cut < 0 || cut + precision + 1 > n_to_s.length() )return n_to_s;
    return n_to_s.substr(0, cut + precision + 1);
}

from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('simulation', views.simulation, name='simulation'),
    # path('about/', views.about, name='about'),
    # path('contact/', views.contact, name='contact'),
]
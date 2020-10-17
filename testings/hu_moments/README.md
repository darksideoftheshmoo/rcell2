Andreas ACL, [23.08.20 13:53]
ahí está. más adelante te puedo ayudar a mejorar el código, pero esto no es el momento para mi 🙂

Andreas ACL, [23.08.20 13:54]
k-means_functions son las funciones (duh!), y k-means_workflow es simplemente un ejemplo de cómo uso el k-means

Andreas ACL, [23.08.20 14:00]
ahora que lograste extraer las máscaras hay una cosa muy buena que podríamos hacer. lo hice antes con las máscaras "feas", pero como no eran perfectas esas máscaras no lo terminé usando

hay unos cálculos morfológicos que se llaman "Hu moments" que me mostró Agustín. son unos cálculos muy locos que juntos terminan representando muy bien la forma de las células. para clasificar, eso es lo más.

yo tengo scripts y todo para hacerlo, solo habría que adaptarlo un poco. podría ser como un post-script luego de importar con R/Cell-ID, que agregaría 8 columnas a la tabla de Rcell, una columna por cada Hu moment. eso después se puede usar para el k-means, o para clasificaciones más copadas aún.

Andreas ACL, [23.08.20 14:41]
es genial. es un poco robando lo que hizo Agustín, pero sinceramente su implementación para clasificar no la íbamos a usar nunca. Además ya hablé con él que le iba a usar su idea y no tuvo ningún problema =)